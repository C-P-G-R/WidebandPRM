import os
import subprocess
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QApplication, QWidget, QTableWidgetItem, QTableWidget, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QScrollArea, QMainWindow
from PyQt5 import QtWidgets, QtCore, QtGui, uic, Qt
from PyQt5.Qt import *
import math
import random
import numpy as np
import pandas as pd
from numpy import trapz
from PIL import Image
import xlsxwriter
from lmfit import models
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from PyQt5.QtGui import QFont
import xml.etree.ElementTree as ET
from scipy.stats import maxwell
import base64
import struct
from collections import defaultdict
from typing import List, Tuple
import csv

class MzMLSpectrum:
    def __init__(self):
        self.scan_num = None
        self.ms_level = None
        self.base_peak_mz = None
        self.base_peak_intensity = None
        self.rt = None
        self.filter_string = None
        self.injection_time = None
        self.isolation_target_mz = None
        self.isolation_lower_offset = None
        self.isolation_upper_offset = None
        self.selected_ion_mz = None
        self.selected_ion_cs = None
        self.selected_ion_intensity = None
        self.scan_description = None
        self.trailer_mz = None
        self.spectrum_ref_scan_num = None
        self.mz_list = None
        self.intensity_list = None
        
class MzMLParser:
    def __init__(self, mzml_file_path):
        self.mzml_file_path = mzml_file_path
        self.spectrum_list = []
        self.spectrum_table = {}
        self.scan_rt_table = {}

    def parse(self):
        print('parse mzML ....')        
        
        with open(self.mzml_file_path, "r") as f:
            xml_data = f.read()

        root = ET.fromstring(xml_data)
        
        namespace = {'ms' : 'http://psi.hupo.org/ms/mzml'}

        spectrum_list = []
        for spectrum_elem in root.findall(".//ms:spectrum", namespace):
            
            spectrum = MzMLSpectrum()
            spectrum.scan_num = spectrum_elem.get("id").split("scan=")[1]
            for param_elem in spectrum_elem.findall(".//ms:cvParam", namespace):
                accession = param_elem.get("accession")
                value = param_elem.get("value")
                if accession == "MS:1000511":
                    spectrum.ms_level = value
                elif accession == "MS:1000504":
                    spectrum.base_peak_mz = value
                elif accession == "MS:1000505":
                    spectrum.base_peak_intensity = value
                elif accession == "MS:1000016":
                    spectrum.rt = float(value)
                elif accession == "MS:1000512":
                    spectrum.filter_string = value
                    spectrum.filter_mz = self.parse_filter_string(spectrum.filter_string)                    
                elif accession == "MS:1000927":
                    spectrum.injection_time = value
                elif accession == "MS:1000827":
                    spectrum.isolation_target_mz = value
                elif accession == "MS:1000828":
                    spectrum.isolation_lower_offset = value
                elif accession == "MS:1000829":
                    spectrum.isolation_upper_offset = value
                elif accession == "MS:1000744":
                    spectrum.selected_ion_mz = value
                elif accession == "MS:1000041":
                    spectrum.selected_ion_cs = value
                elif accession == "MS:1000042":
                    spectrum.selected_ion_intensity = value
            
            for user_param_elem in spectrum_elem.findall(".//ms:userParam", namespace):
                name = user_param_elem.get("name")
                value = user_param_elem.get("value")
                if name == "scan description":
                    spectrum.scan_description = value
                    spectrum.scan_description_cs = self.extract_scan_description_cs(spectrum.scan_description)
                    if spectrum.scan_description.find("_") != -1:
                        spectrum.spectrum_type = 2
                    elif spectrum.ms_level == 2:
                        spectrum.spectrum_type = 1
                    elif spectrum.ms_level == 1:
                        spectrum.spectrum_type = 0
                elif name == "[Thermo Trailer Extra]Monoisotopic M/Z:":
                    spectrum.trailer_mz = value

            binary_elems = spectrum_elem.findall(".//ms:binary", namespace)
            
            for i, binary in enumerate(binary_elems):
                encoded_bytes = binary.text
                peaks = base64.b64decode(encoded_bytes)
                byte_order = '<'  # Little-endian
                
                if i == 0:
                    mz_values = struct.unpack(byte_order + 'f' * (len(peaks) // 4), peaks)
                    spectrum.mz_list = [round(mz, 4) for mz in mz_values]
                else:
                    intensity_values = struct.unpack(byte_order + 'f' * (len(peaks) // 4), peaks)
                    spectrum.intensity_list = [round(intensity) for intensity in intensity_values]
            
            self.spectrum_list.append(spectrum)
            self.spectrum_table[spectrum.scan_num] = spectrum
            self.scan_rt_table[spectrum.scan_num] = spectrum.rt
        
        return self.spectrum_list, self.spectrum_table, self.scan_rt_table   
    
    def parse_filter_string(self, filter_string):
        # Find the index of the first occurrence of the string "ms2"
        ms2_index = filter_string.find("ms2")
    
        # Find the index of the first occurrence of the character "@"
        at_index = filter_string.find("@")
    
        # Extract the substring between the two indices
        mz_str = filter_string[ms2_index + 4:at_index]
    
        # Convert the substring to a float
        filter_mz = float(mz_str)
    
        return filter_mz 
    def extract_scan_description_cs(self, scan_description):
        # Extract the last character from the scan description string
        cs_str = scan_description[-1]
        # Convert the last character to an integer
        scan_description_cs = int(cs_str)
        # Return the extracted scan description CS
        return scan_description_cs

class PRMTarget:
    def __init__(self):
        self.peptide = ""
        self.mz = ""
        self.rt = ""
        self.window = ""
        self.group_id = ""
        self.transitions = list()
        self.transitions_endo = list()
        self.all_transition_profile = list()
        self.all_transition_profile_endo = list()

class PRMTransition:
    def __init__(self):
        self.ion_type = ""
        self.mz = ""
        self.group_id = ""
        self.profile = []

class WidebandPRM:
    def __init__(self, args, spectrum_list, spectrum_table, scan_rt_table):
        self.spectrum_list = spectrum_list
        self.spectrum_table = spectrum_table
        self.scan_rt_table = scan_rt_table
        
        self.target1_url = ""
        self.target2_url = ""
        self.target3_url = ""
        self.trigger1_url = ""
        self.trigger2_url = ""
        self.trigger3_url = ""
        self.group_id_list = []
        self.mzml_url = ""
        self.precursor_tolerance = 10.0
        self.msms_tolerance = 10.0
        self.peak_count = 100
        self.is_logging = False
        self.target_table = defaultdict(list)
        self.logging_target_peptide = "a"

        self.heavy_form = 6.02013  # +1
        self.heavy_form_2 = 3.01006  # +2
        self.heavy_form_3 = 2.00671  # + 3
        self.heavy_form_4 = 1.50503  # + 4

        self.parse_args(args)

    def parse_args(self, args):
        i = 0
        while i < len(args):
            if args[i].upper() == "-TARGET1":
                self.target1_url = args[i + 1]
                self.trigger1_url = os.path.splitext(self.target1_url)[0] + "_Trigger.csv"
            elif args[i].upper() == "-TARGET2":
                self.target2_url = args[i + 1]
                self.trigger2_url = os.path.splitext(self.target2_url)[0] + "_Trigger.csv"
            elif args[i].upper() == "-MZML":
                self.mzml_url = args[i + 1]
            elif args[i].upper() == "-PEAKCOUNT":
                self.peak_count = int(args[i + 1])
            elif args[i].upper() == "-TOLERANCE":
                self.precursor_tolerance = float(args[i + 1])
                self.msms_tolerance = float(args[i + 1])
            elif args[i].upper() == "-V":
                self.is_logging = True
            i += 2

    def execute(self):
        print("##### Start #####")
        print(f"Precursor m/z Tolerance : {self.precursor_tolerance}")
        print(f"MS2 m/z Tolerance : {self.msms_tolerance}")
        print(f"No. Peaks for Quantify : {self.peak_count}")

        self.load_target(self.target1_url)
        self.load_trigger(self.trigger1_url)

        if self.target2_url and self.target2_url != "":
            self.load_target(self.target2_url)
            self.load_trigger(self.trigger2_url)

        for spectrum in self.spectrum_list:
            
            if spectrum.spectrum_type == 2:  # Quant-Mode spectrum
                target = self.get_target(spectrum)
                if target is None:
                    continue
                
                sum_intensity = self.get_sum_of_transition(spectrum, target, spectrum.mz_list, spectrum.intensity_list)

                rt_intensity = [spectrum.rt, sum_intensity[0], spectrum.scan_num]
                target.all_transition_profile.append(rt_intensity)
                
                rt_intensity_endo = [spectrum.rt, sum_intensity[1], spectrum.scan_num]
                target.all_transition_profile_endo.append(rt_intensity_endo)

        output = self.write_result()        
        print("mzML loaded!")
        return output

    def get_target(self, spectrum):
        filter_mz = spectrum.filter_mz
        scan_description_cs = spectrum.scan_description_cs

        find_mz = 0.0

        if scan_description_cs == 2:
            find_mz = filter_mz + 1.5
        elif scan_description_cs == 3:
            find_mz = filter_mz + 1.0
        elif scan_description_cs == 4:
            find_mz = filter_mz + 0.75
        else:
            raise Exception(f"Unknown Scan Description CS {spectrum.scan_num}")

        target_candidates = []

        for group_id in self.group_id_list:
            target_list = self.target_table[group_id]
            for target in target_list:
                mma = self.get_ppm(find_mz, target.mz)

                if mma < self.precursor_tolerance:
                    rt_from = target.rt - target.rt_tolerance
                    rt_to = target.rt + target.rt_tolerance

                    if rt_from <= spectrum.rt <= rt_to:
                        target_candidates.append(target)


        if len(target_candidates) == 0:
            return None

        if len(target_candidates) == 1:
            return target_candidates[0]

        # Get the nearest target from candidates
        highest_target = None
        most_intensity = 0
        # 후보가 여러개라면 transition과 매칭되는 개수를 구해서 가장 많이 매칭되는 것을 선택        
        for target in target_candidates:
            transitions = target.transitions
            for tr in transitions:
                intensity = self.get_peak_intensity(spectrum, tr, spectrum.mz_list, spectrum.intensity_list)
                if most_intensity < intensity:
                    most_intensity = intensity
                    highest_target = target

        return highest_target

    def load_target(self, target_url):
        print(f"Load target : {target_url}")

        with open(target_url, newline="", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            target = PRMTarget()

            for row in reader:
                target = PRMTarget()
                target.peptide = row[0]
                target.mz = float(row[1])
                target.rt = float(row[3])
                target.window = float(row[4])
                target.rt_tolerance = target.window / 2.0
                group_id = target.mz

                target_list = self.target_table[group_id]          
                target_list.append(target)
                self.group_id_list.append(group_id)
            print("No. target : ", len(self.target_table))

    def load_trigger(self, trigger_url):
        print(f"Load trigger : {trigger_url}")

        with open(trigger_url, newline="", encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            trigger = PRMTransition()

            transitions = list()
            previousGroupId = ''
            transition = None
            endo_transition = None
            target = None
            
            
            for row in reader:
                
                if previousGroupId != '' and previousGroupId != row[2]:
                    target_list = self.target_table[float(previousGroupId)]
                    for i in range(len(target_list)):
                        target = target_list[i]
                        if len(target.transitions) == 0:
                            target.transitions = transitions
                            
                            # endo list append
                            for tr in target.transitions:
                                endo_transition = PRMTransition()
                                endo_mz = self.calculate_endo_mz(tr.mz, tr.ion_type)
                                endo_ion_type = tr.ion_type
                                endo_transition.mz = endo_mz
                                endo_transition.group_id = tr.group_id
                                target.transitions_endo.append(endo_transition)
                            transitions = list()
                            break
                        
                transition = PRMTransition()
                transition.ion_type = row[0]
                transition.mz = float(row[1])
                transition.group_id = row[2]
                transitions.append(transition)
                
                previousGroupId = row[2]
                            
            target_list = self.target_table[float(previousGroupId)]
            for i in range(len(target_list)):
                target = target_list[i]
                if len(target.transitions) == 0:
                    target.transitions = transitions
                    # endo list append
                    for tr in target.transitions:
                        endo_transition = PRMTransition()
                        endo_mz = self.calculate_endo_mz(tr.mz, tr.ion_type)
                        endo_ion_type = tr.ion_type
                        endo_transition.mz = endo_mz
                        endo_transition.group_id = tr.group_id
                        target.transitions_endo.append(endo_transition)
                    transitions = list()
                    break
            
    def calculate_endo_mz(self, tr_mz, ion_type):
        
        HEAVY_FORM = 6.02013
        HEAVY_FORM_2 = 3.010065
        HEAVY_FORM_3 = 2.00671
        HEAVY_FORM_4 = 1.50503
        
        endo_mz = 0.0
        cs_str = ion_type[-1]
        cs = int(cs_str)
    
        if cs == 1:
            endo_mz = tr_mz - HEAVY_FORM
        elif cs == 2:
            endo_mz = tr_mz - HEAVY_FORM_2
        elif cs == 3:
            endo_mz = tr_mz - HEAVY_FORM_3
        elif cs == 4:
            endo_mz = tr_mz - HEAVY_FORM_4
        else:
            endo_mz = tr_mz - HEAVY_FORM_2

        return endo_mz

    def get_sum_of_transition(self, spectrum, target, mz_list, intensity_list):
        sum_intensity = 0
        sum_intensity_endo = 0
        
        transition_list = target.transitions
        transition_list_endo = target.transitions_endo
        
        for i in range(0, len(transition_list)):
            # SIL
            transition = transition_list[i]
                        
            if not transition.ion_type.startswith("y"):
                continue
            
            intensity = self.get_peak_intensity(spectrum, transition, mz_list, intensity_list)
            sum_intensity += intensity
            
            # ENDO
            transition_endo = transition_list_endo[i]
            intensity_endo = self.get_peak_intensity(spectrum, transition_endo, mz_list, intensity_list)
            sum_intensity_endo += intensity_endo
            
            profile = [str(spectrum.scan_num), str(spectrum.rt), str(intensity)]
            transition.profile.append(profile)  # transition 각각의 profile을 저장 
            
            profile_endo = [str(spectrum.scan_num), str(spectrum.rt), str(intensity_endo)]
            transition_endo.profile.append(profile_endo)  # transition 각각의 profile을 저장 
            
            if (i + 1) >= self.peak_count:                
                break
        
        return [sum_intensity, sum_intensity_endo]
    
    
    def get_peak_intensity(self, spectrum, transition, mz_list, intensity_list):
        intensity = 0
        
        index = int(len(mz_list) / 2)
        index_array = list()
        index = self.find_peak(spectrum, transition.mz, mz_list, index, index_array)
        
        found_mz = mz_list[index]
        intensity_of_index = intensity_list[index]
        
        most_min_diff = self.get_ppm(mz_list[index], transition.mz)
        
        most_min_diff_mz = 9999999.0
        most_min_diff_intensity = 0
        most_min_diff_peak_index = -1
        
        # 앞으로 이동 하면서 3ppm 이내인 peak들 수집
        for i in range(index - 1, 0, -1):
            mma = self.get_ppm(mz_list[i], transition.mz)
            if mma > 10.0:
                break
            else:
                if most_min_diff > mma:
                    most_min_diff = mma
                    most_min_diff_intensity = intensity_list[i]
                    most_min_diff_peak_index = i
        
        # 뒤로 이동 하면서 3ppm 이내인 peak들 수집
        for i in range(index + 1, len(mz_list)):
            mma = self.get_ppm(mz_list[i], transition.mz)
            if mma > 10.0:
                break
            else:
                if most_min_diff > mma:
                    most_min_diff = mma
                    most_min_diff_intensity = intensity_list[i]
                    most_min_diff_peak_index = i
        
        # mostAbundantPeakIndexrk -1이면 근처에 peak이 없는 경우이니깐 index 위치를 리턴
        # 아니면 abundance를 비교해서 높은값을 가지는 index를 리턴
        if most_min_diff_peak_index == -1:
            mma = self.get_ppm(mz_list[index], transition.mz)
            if mma > 10.0:
                return 0
            else:                
                return intensity_of_index
        else:
            mma = self.get_ppm(mz_list[most_min_diff_peak_index], transition.mz)
            if mma > 10.0:
                return 0
            else:
                return most_min_diff_intensity
    
    
    def find_peak(self, spectrum, find_mz, mz_list, index, index_array):
        mz = mz_list[index]
        
        mma = self.get_ppm(mz, find_mz)
        if mma <= 10.0:
            return index
        
        if mz > find_mz:
            new_index = index - int((index - self.find_lower_val(index_array, index)) / 2)
            if index - new_index == 0:
                return index
            index_array.append(index)
            return self.find_peak(spectrum, find_mz, mz_list, new_index, index_array)
        else:
            new_index = index + int((self.find_upper_val(index_array, index, len(mz_list)) - index) / 2)
            if new_index - index == 0:
                return index
            index_array.append(index)
            return self.find_peak(spectrum, find_mz, mz_list, new_index, index_array)
    
    def find_upper_val(self, index_array, index, limit_index):
        upper_index = -1
        for i in range(len(index_array) - 1, -1, -1):
            pre_index = index_array[i]
            if pre_index > index:
                upper_index = pre_index
                break
        if upper_index == -1:
            upper_index = limit_index
        return upper_index
    
    def find_lower_val(self, index_array, index):
        lower_index = -1
        for i in range(len(index_array) - 1, -1, -1):
            pre_index = index_array[i]
            if pre_index < index:
                lower_index = pre_index
                break
        if lower_index == -1:
            lower_index = 0
        return lower_index
    
    def get_ppm(self, mz1, mz2):
        return abs((mz1 - mz2) / mz2) * 1e6

    def write_result(self):
        
        output = self.mzml_url[0:self.mzml_url.rfind(".")] + "_PeakTracing.txt"
        print(f"Write result : {output}")

        with open(f"{output}", mode="w", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')

            header = ["Peptide","Precursor(m/z)","ScanNum","IonType","Transition(m/z)","RT","Intensity(SIL)","Intensity(Endo)","Ratio(Endo/SIL)"]

            writer.writerow(header)

            for group_id in self.group_id_list:
                target_list = self.target_table[group_id]

                for target in target_list:
                    
                    transition_list = target.transitions
                    transition_list_endo = target.transitions_endo
                    
                    for i in range (0, len(transition_list)):
                        
                        transition = transition_list[i]
                        transition_endo = transition_list_endo[i]
                        
                        if transition.ion_type[0:1] != 'y':
                            continue
                        
                        profile_sil = transition.profile
                        profile_endo = transition_endo.profile
                        
                        for j in range (0, len(profile_sil)):
                            point_sil = profile_sil[j]
                            point_endo = profile_endo[j]
                            
                            ratio = 'NaN'
                            if float(point_sil[2]) > 0:
                                ratio = float(point_endo[2]) / float(point_sil[2])
                                ratio = round(ratio, 4)
                            
                            row = list()
                            row.append(target.peptide)
                            row.append(target.mz)
                            row.append(point_sil[0])
                            row.append(transition.ion_type)
                            row.append(transition.mz)
                            row.append(point_sil[1])
                            row.append(point_sil[2])
                            row.append(point_endo[2])
                            row.append(ratio)
                            
                            writer.writerow(row)
        return output

class Spectrum:
    def __init__(self):
        self.scan_num = 0
        self.spectrum_type = 0
        self.scan_description_cs = 0
        self.filter_string_mz = 0.0
        self.rt = 0.0
        self.mz_list = []
        self.intensity_list = []

class BoltzmannModel:
    
    def __init__(self, model_x, model_y, model_sigma, model_center):
        self.model_x = model_x
        self.model_y = model_y
        self.model_sigma = model_sigma
        self.model_center = model_center
        
    def get_model_sigma(self):
        return self.model_sigma
        
    def get_x(self):
        return self.model_x
    
    def get_y(self):
        return self.model_y
    
    def get_model_center(self):
        return self.model_center
    
    def get_model_area(self):        
        return trapz(self.model_y, self.model_x)
    
    def get_model_area_with_range(self, min_x, max_x):
        #print('get_model_area_with_range():', min_x, '~', max_x)
        new_model_y = list()
        new_model_x = list()
        for i in range (len(self.model_x)):
            if self.model_x[i] >= min_x and self.model_x[i] <= max_x:
                new_model_x.append(self.model_x[i])
                new_model_y.append(self.model_y[i])
        return trapz(new_model_y, new_model_x)
    
    def set_sse(self, sse):
        self.sse = sse
    
    def get_sse(self):
        return self.sse
    
    def get_peak_width(self):
        min_x = min(self.model_x)
        max_x = max(self.model_x)
        peak_width = max_x - min_x
        peak_width = round(peak_width, 2)
        return peak_width
    
    def get_model_range(self):
        #print('call get_model_range() ...')
        min_model_x = min(self.model_x)
        max_model_x = max(self.model_x)
        #print("min-max:", min_model_x, "~", max_model_x)
        return min_model_x, max_model_x
    
class Transition:
        
    def __init__(self, precursorMz, ionType):
        self.precursorMz = precursorMz
        self.ionType = ionType
        self.profileSIL = {}
        self.profileEndo = {}
        self.model_SIL = list()
        self.model_Endo = list()        
        self.peak_width_list = list()
        self.model_cnt_SIL = 1
        self.model_cnt_Endo = 1
        self.final_model_index_SIL = 0
        self.final_model_index_Endo = 0
        self.rt_min_for_view = -1
        self.rt_max_for_view = -1
        self.use_transition = 0 #0:ratio 구할 때 미사용  1:사용 
        self.used_replicate_1 = 0 # replicate #1에서 사용했었는지 
        self.used_replicate_2 = 0 # replicate #2에서 사용했었는지 
        self.rawMinRT = 999999
        self.rawMaxRT = 0
        
    def set_img_path(self, img_path):
        self.img_path = img_path
    def get_img_path(self):
        return self.img_path
        
    def set_used_replicate_1(self, replicate_1):
        self.used_replicate_1 = replicate_1
    def get_used_replicate_1(self):
        return self.used_replicate_1
    
    def set_used_replicate_2(self, replicate_2):
        self.used_replicate_2 = replicate_2
    def get_used_replicate_2(self):
        return self.used_replicate_2
        
    def set_use_transition(self, use_transition):
        self.use_transition = use_transition
    def get_use_transition(self):
        return self.use_transition
        
    def set_rt_min_for_view(self, rt_min):
        self.rt_min_for_view = rt_min
    def get_rt_min_for_view(self):
        return self.rt_min_for_view
    
    def set_rt_max_for_view(self, rt_max):
        self.rt_max_for_view = rt_max
    def get_rt_max_for_view(self):
        return self.rt_max_for_view
        
    def get_model_SIL(self):
        return self.model_SIL
    
    def get_model_Endo(self):
        return eslf.model_Endo
        
    def set_final_model_index_SIL(self, index):
        self.final_model_index_SIL = index
    
    def set_final_model_index_Endo(self, index):
        self.final_model_index_Endo = index
        
    def get_final_model_index_SIL(self):
        return self.final_model_index_SIL
    
    def get_final_model_index_Endo(self):
        return self.final_model_index_Endo
        
    def add_model_SIL(self, model):
        index = -1
        for i in range(0, len(self.model_SIL)):
            eModel = self.model_SIL[i]
            if eModel.get_model_center() > model.get_model_center():
                index = i
                break
        
        if index == -1:
            self.model_SIL.append(model)            
            #print("adding model = -1 :", model.get_model_center())
        else:
            self.model_SIL.insert(index, model)
            #print("adding model != -1 :", model.get_model_center())
    
    def add_model_Endo(self, model):
        index = -1
        for i in range(0, len(self.model_Endo)):
            eModel = self.model_Endo[i]
            if eModel.get_model_center() > model.get_model_center():
                index = i
                break
        
        if index == -1:
            self.model_Endo.append(model)
        else:
            self.model_Endo.insert(index, model)
            
    def get_peak_width(self):
        return max(self.peak_width_list)
    
    def getPrecursorMz(self):
        return self.precursorMz
    def setPrecursorMz(self, precursorMz):
        self.precursorMz= precursorMz
        
    def setTransitionMz(self, transitionMz):
        self.transitionMz = transitionMz
    def getTransitionMz(self):
        return self.transitionMz
        
    def setIonType(self, ionType):
        self.ionType = ionType
    def getIonType(self):
        return self.ionType
    
    def setProfileSIL(self, profileSIL):
        self.profileSIL = profileSIL
    
    def setProfileEndo(self, profileEndo):
        self.profileEndo = profileEndo
    
    def getProfileSIL(self):
        return self.profileSIL
    
    def getProfileEndo(self):
        return self.profileEndo
    
    def getProfileSIL(self):
        return self.profileSIL
    
    def getProfileEndo(self):
        return self.profileEndo
    
    def setScanNum(self, scanNum):
        self.scanNum = scanNum
    def getScanNum(self):
        return self.scanNum
    
    def setRt(self, rt):
        self.rt = rt
    def getRt(self):
        return self.rt
    
    def setImgPath(self, imgPath):
        self.imgPath = imgPath
    def getImgPath(self):
        return self.imgPath
    
    def setRawMinRT(self, minRt):
        self.rawMinRT = minRt
    def getRawMinRT(self):
        return self.rawMinRT
    
    def setRawMaxRT(self, maxRt):
        self.rawMaxRT = maxRt
    def getRawMaxRT(self):
        return self.rawMaxRT
    
    # endo가 SIL보다 10배 이상 크면 1을 리턴함으로써 Endo 우선 modeling 구현 
    def checkEndoPriority(self):
        intensityListSIL = list(self.profileSIL.values())
        intensityListEndo = list(self.profileEndo.values())
        sil_max_intensity = max(intensityListSIL)
        endo_max_intensity = max(intensityListEndo)
        
        if sil_max_intensity * 5 < endo_max_intensity:
            return 1
        else:
            return 0       
    
    def fit_model(self):
        rtList = list(self.profileSIL.keys())
        self.setRawMinRT(min(rtList))
        self.setRawMaxRT(max(rtList))
        intensityListSIL = list(self.profileSIL.values())
        intensityListEndo = list(self.profileEndo.values())
        
        #print(self.getIonType(), ":", len(rtList), ":", len(intensityListSIL), ",", len(intensityListEndo))
        
        if len(rtList) < 3:
            return -1
        
        new_x_SIL, new_y_SIL, sigma, center = self.make_model_SIL(rtList, intensityListSIL)
        new_x_Endo, new_y_Endo = self.make_model_Endo(rtList, intensityListEndo, sigma, center)
        
        return 1

    def fit_model_reverse(self):
        rtList = list(self.profileSIL.keys())
        self.setRawMinRT(min(rtList))
        self.setRawMaxRT(max(rtList))
        intensityListSIL = list(self.profileSIL.values())
        intensityListEndo = list(self.profileEndo.values())
        
        if len(rtList) < 5:
            return -1
        
        new_x_Endo, new_y_Endo, sigma, center = self.make_model_Endo_Reverse(rtList, intensityListEndo)
        new_x_SIL, new_y_SIL = self.make_model_SIL_Reverse(rtList, intensityListSIL, sigma, center)
        
        return 1
    
    def fit_model_custom(self, model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, m_factor):
        
        self.model_SIL = list()
        self.model_Endo = list()
        
        self.model_cnt_SIL = model_cnt_SIL
        self.model_cnt_Endo = model_cnt_Endo
        self.rt_from = rt_from
        self.rt_to = rt_to
        self.m_factor = m_factor
        
        #print(self.ionType, " modeling ...")
        rtListORG = list(self.profileSIL.keys())
        intensityListSIL = list(self.profileSIL.values())
        rtList, intensityListSIL = self.check_range(rtListORG, intensityListSIL, rt_from, rt_to)
        
        intensityListEndo = list(self.profileEndo.values())        
        rtList, intensityListEndo = self.check_range(rtListORG, intensityListEndo, rt_from, rt_to)
        
        #print("fit_model_custom:SIL", len(rtList), ",", len(intensityListSIL), " Endo-", len(rtList), ",", len(intensityListEndo))
        
        if len(rtList) < 5:
            return -1
        
        new_x_SIL, new_y_SIL, sigma, center = self.make_model_SIL(rtList, intensityListSIL)
        new_x_Endo, new_y_Endo = self.make_model_Endo(rtList, intensityListEndo, sigma, center)
        return 1
        
    def fit_model_custom_reverse(self, model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, m_factor):
        
        self.model_SIL = list()
        self.model_Endo = list()
        
        self.model_cnt_SIL = model_cnt_SIL
        self.model_cnt_Endo = model_cnt_Endo
        self.rt_from = rt_from
        self.rt_to = rt_to
        self.m_factor = m_factor
        
        #print(self.ionType, " modeling ...")
        rtListORG = list(self.profileEndo.keys())
        
        intensityListEndo = list(self.profileEndo.values())
        rtList, intensityListEndo = self.check_range(rtListORG, intensityListEndo, rt_from, rt_to)
        
        intensityListSIL = list(self.profileSIL.values())        
        rtList, intensityListSIL = self.check_range(rtListORG, intensityListSIL, rt_from, rt_to)
        
        if len(rtList) < 5:
            return -1
        
        new_x_Endo, new_y_Endo, sigma, center = self.make_model_Endo_Reverse(rtList, intensityListEndo)
        new_x_SIL, new_y_SIL = self.make_model_SIL_Reverse(rtList, intensityListSIL, sigma, center)
        return 1       

    # 해당 메소드는 사용자 입력을 받는 부분이므로 더 이상 사용하지 않음         
    def set_contact(self):
        is_retry = 1
        while True:
            is_retry = input("retry modeling?")
            if is_retry.strip() == "0" or is_retry.strip() == "1":
                break
            continue
            
        return int(is_retry)

    def check_range(self, rtList, intensityList, min_x, max_x):
        new_rt_list = list()
        new_intensity_list = list()
        
        for i in range(0, len(rtList)):
            rt = rtList[i]
            if rt < min_x or rt > max_x:
                continue
            new_rt_list.append(rt)
            new_intensity_list.append(intensityList[i])
        
        return new_rt_list, new_intensity_list
    
    def getMovingAverage(self, rtList, intensityList, w):
        numpyArray = np.array(intensityList)
        m_avg_values = np.convolve(numpyArray, np.ones(w), 'valid') / w
        m_avg_rt_intensity = {}
        for i in range(0, len(rtList)):
            rt = rtList[i]
            intensity = 0
            if i < (w-1):
                intensity = intensityList[i]
            else: 
                intensity = m_avg_values[i - (w-1)]
            
            m_avg_rt_intensity[rt] = intensity
        return m_avg_rt_intensity 
        
    def make_model_SIL(self, rt, intensity):
        
        self.model_SIL = list()
        #MODEL_TYPE = "GaussianModel"
        MODEL_TYPE = "LorentzianModel"
        #MODEL_TYPE = "VoigtModel"
        spec = {
            'x': rt,
            'y': intensity,
            'model': [
                {'type': MODEL_TYPE} for i in range(0, self.model_cnt_SIL)
            ]
        }
        model, params = self.generate_model(spec)
        output = model.fit(spec['y'], params, x=spec['x'])
        
        for i in range(0, self.model_cnt_SIL):
            model_index= "m" + str(i) + "_"
            center = output.params[model_index + 'center'].value
            sigma = output.params[model_index + 'sigma'].value
            height = output.params[model_index + 'height'].value
            #amp    = output.params[model_index + 'amplitude'].value
            #fwhm   = output.params[model_index + 'fwhm'].value
            
            new_x, new_y, sse = self.doBoltzmann(rt, intensity, center, sigma, height, 0)
            #new_x, new_y, sse = self.doBoltzmann_refine(rt, intensity, center, sigma, height, 0)
            
            boltzmannModel = BoltzmannModel(new_x, new_y, sigma, center)
            boltzmannModel.set_sse(sse)
            self.add_model_SIL(boltzmannModel)
            
        return new_x, new_y, sigma, center
    
    def make_model_Endo(self, rt, intensity, sigmaSIL, centerSIL):
        self.model_Endo = list()
        #MODEL_TYPE = "GaussianModel"
        MODEL_TYPE = "LorentzianModel"
        #MODEL_TYPE = "VoigtModel"
        spec = {
            'x': rt,
            'y': intensity,
            'model': [
                {'type': MODEL_TYPE} for i in range(0, self.model_cnt_Endo)
            ]
        }
        model, params = self.generate_model_with_hint(spec,sigmaSIL, centerSIL)
        output = model.fit(spec['y'], params, x=spec['x'])
        
        for i in range(0, self.model_cnt_Endo):
            model_index= "m" + str(i) + "_"
            center = output.params[model_index + 'center'].value
            sigma = output.params[model_index + 'sigma'].value
            height = output.params[model_index + 'height'].value
            
            new_x, new_y, sse = self.doBoltzmann(rt, intensity, center, sigma, height, 0)
            #new_x, new_y, sse = self.doBoltzmann_refine(rt, intensity, center, sigma, height, 0)
            boltzmannModel = BoltzmannModel(new_x, new_y, sigma, center)
            boltzmannModel.set_sse(sse)
            self.add_model_Endo(boltzmannModel)
            
        return new_x, new_y
    
    def make_model_SIL_Reverse(self, rt, intensity, sigmaSIL, centerSIL):
        self.model_SIL = list()
        #MODEL_TYPE = "GaussianModel"
        MODEL_TYPE = "LorentzianModel"
        #MODEL_TYPE = "VoigtModel"
        spec = {
            'x': rt,
            'y': intensity,
            'model': [
                {'type': MODEL_TYPE} for i in range(0, self.model_cnt_SIL)
            ]
        }
        model, params = self.generate_model_with_hint(spec,sigmaSIL, centerSIL)
        output = model.fit(spec['y'], params, x=spec['x'])
        
        for i in range(0, self.model_cnt_SIL):
            model_index= "m" + str(i) + "_"
            center = output.params[model_index + 'center'].value
            sigma = output.params[model_index + 'sigma'].value
            height = output.params[model_index + 'height'].value
            
            new_x, new_y, sse = self.doBoltzmann(rt, intensity, center, sigma, height, 0)
            #new_x, new_y, sse = self.doBoltzmann_refine(rt, intensity, center, sigma, height, 0)
            boltzmannModel = BoltzmannModel(new_x, new_y, sigma, center)
            boltzmannModel.set_sse(sse)
            self.add_model_SIL(boltzmannModel)
            #print("add SIL Model ", self.ionType)
        return new_x, new_y
    
    def make_model_Endo_Reverse(self, rt, intensity):
        
        self.model_Endo = list()
        #MODEL_TYPE = "GaussianModel"
        MODEL_TYPE = "LorentzianModel"
        #MODEL_TYPE = "VoigtModel"
        spec = {
            'x': rt,
            'y': intensity,
            'model': [
                {'type': MODEL_TYPE} for i in range(0, self.model_cnt_Endo)
            ]
        }
        model, params = self.generate_model(spec)
        output = model.fit(spec['y'], params, x=spec['x'])
        
        for i in range(0, self.model_cnt_Endo):
            model_index= "m" + str(i) + "_"
            center = output.params[model_index + 'center'].value
            sigma = output.params[model_index + 'sigma'].value
            height = output.params[model_index + 'height'].value
            #amp    = output.params[model_index + 'amplitude'].value
            #fwhm   = output.params[model_index + 'fwhm'].value
            
            new_x, new_y, sse = self.doBoltzmann(rt, intensity, center, sigma, height, 0)
            #new_x, new_y, sse = self.doBoltzmann_refine(rt, intensity, center, sigma, height, 0)
            
            boltzmannModel = BoltzmannModel(new_x, new_y, sigma, center)
            boltzmannModel.set_sse(sse)
            self.add_model_Endo(boltzmannModel)
            #print("add Edno Model ", self.ionType)
        return new_x, new_y, sigma, center
    ######################## 모델 생성 ##########################
    def generate_model(self, spec):
        
        composite_model = None
        params = None
        
        x = spec['x']
        y = spec['y']
        x_min = np.min(x)
        x_max = np.max(x)
        x_range = x_max - x_min
        y_max = np.max(y)
        
        for i, basis_func in enumerate(spec['model']):
            
            prefix = f'm{i}_'
            model = getattr(models, basis_func['type'])(prefix=prefix)
            if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
                model.set_param_hint('sigma', min=1e-6, max=x_range)
                #model.set_param_hint('sigma', min=1e-3, max=x_range)
                model.set_param_hint('center', min=x_min, max=x_max)                
                model.set_param_hint('height', min=1e-6, max=1.1*y_max)
                #model.set_param_hint('height', min=1e-3, max=1.1*y_max)
                model.set_param_hint('amplitude', min=1e-6)
                #model.set_param_hint('amplitude', min=1e-3)
                # default guess is horrible!! do not use guess()
                default_params = {
                    prefix+'center': x_min + x_range * random.random(),
                    prefix+'height': y_max * random.random(),
                    prefix+'sigma': x_range * random.random()
                }
            else:
                raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
            if 'help' in basis_func:  # allow override of settings in parameter
                for param, options in basis_func['help'].items():
                    model.set_param_hint(param, **options)
            model_params = model.make_params(**default_params, **basis_func.get('params', {}))
            if params is None:
                params = model_params
            else:
                params.update(model_params)
            if composite_model is None:
                composite_model = model
            else:
                composite_model = composite_model + model
        return composite_model, params    
    
    ######################## 모델 생성 ##########################
    def generate_model_with_hint(self, spec, sigma, center):
        
        composite_model = None
        params = None
        
        x = spec['x']
        y = spec['y']
        x_min = np.min(x)
        x_max = np.max(x)
        x_range = x_max - x_min
        y_max = np.max(y)
        sigma_min = sigma - 0.01
        sigma_max = sigma + 0.01
        center_min = center - 0.03
        center_max = center + 0.03
        
        #print("sigma:", sigma, " center:", center)
        
        for i, basis_func in enumerate(spec['model']):
            
            prefix = f'm{i}_'
            model = getattr(models, basis_func['type'])(prefix=prefix)
            if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
                #model.set_param_hint('sigma', min=1e-6, max=x_range)
                model.set_param_hint('sigma', min=sigma_min, max=sigma_max)
                #model.set_param_hint('center', min=x_min, max=x_max)                
                model.set_param_hint('center', min=center_min, max=center_max)                
                model.set_param_hint('height', min=1e-6, max=1.1*y_max)
                model.set_param_hint('amplitude', min=1e-6)

                # default guess is horrible!! do not use guess()
                default_params = {
                    prefix+'center': center,
                    prefix+'height': y_max * random.random(),
                    prefix+'sigma': sigma
                }
                '''
                default_params = {
                    prefix+'center': x_min + x_range * random.random(),
                    prefix+'height': y_max * random.random(),
                    prefix+'sigma': x_range * random.random()
                }
                '''
            else:
                raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
            if 'help' in basis_func:  # allow override of settings in parameter
                for param, options in basis_func['help'].items():
                    model.set_param_hint(param, **options)
            model_params = model.make_params(**default_params, **basis_func.get('params', {}))
            if params is None:
                params = model_params
            else:
                params.update(model_params)
            if composite_model is None:
                composite_model = model
            else:
                composite_model = composite_model + model
        return composite_model, params     
    
    def doBoltzmann(self, data_x, data_y, center, sigma, height, baseline):
        sigma = sigma - 0.0045
        #center =center - 0.008
        samplesize = 10000
        max_intensity = height
        loc = 0 #my guess
        a_value = np.sqrt((sigma**2 * math.pi)/(3*math.pi - 8)) #calculated based on wiki description
        v_2d = maxwell.rvs(loc, a_value, size=samplesize) #array corresponding to 2D proper motion obtained from Hubbs
        mean, var, skew, kurt = maxwell.stats(moments='mvsk')
        maxx = np.linspace(min(v_2d), max(v_2d), samplesize)
        params = maxwell.fit(v_2d, floc=0)
        factor_y = max_intensity / np.max(maxwell.pdf(maxx, *params))
        new_y = maxwell.pdf(maxx, *params) * factor_y
        max_var = max(new_y)
        max_idx = -1
        for i in range(0, len(new_y)):
            if new_y[i] == max_var:
                max_idx = i
                break
        max_var_x = maxx[max_idx]
        factor_x = center - max_var_x
        new_x = maxx + factor_x
        
        # SSE 계산 
        shift_rt_arr = list()
        shift_rt_arr.append(0.0)
        min_sse = -1
        min_sse_rt = -1
        for idx in range(0, len(shift_rt_arr)):
            shift_rt = shift_rt_arr[idx]
            
            shifted_rt_arr = self.rt_shift(data_x, shift_rt)
            
            sse_list = list()
            
            for i in range(0, len(shifted_rt_arr)):
                x = shifted_rt_arr[i] - factor_x
                
                x_min = min(new_x)
                x_max = max(new_x)
                if shifted_rt_arr[i] < x_min or x > x_max:
                    continue
                
                _var = maxwell.pdf(x, *params)
                
                y = data_y[i]
                _var_weight = _var * factor_y
                sse = (y - _var_weight) **2
                sse_list.append(sse)
                #print("->\t", x, "\t", sse)

            if idx == 0:
                min_sse = sum(sse_list)
                min_sse_rt = shift_rt
                #print("SSE:", sum(sse_list))
            else:
                if min_sse > sum(sse_list):
                    min_sse = sum(sse_list)
                    min_sse_rt = shift_rt
                    #print("SSE:", sum(sse_list))
                    
        #print("Minimum SSE : ", min_sse_rt)        
        optimized_x = self.rt_shift(new_x, -min_sse_rt)
        
        # 범위 내의 +- 0.3 구간출력 
        refined_x = list();
        refined_y = list();
        
        #min_x = np.min(data_x) - 0.3
        #max_x = np.max(data_x) + 0.3
        min_x = np.min(new_x) - 0.1
        max_x = np.max(new_x) + 0.1
        
        for i in range(0, len(new_x)):
            rt = new_x[i]
            intensity = new_y[i]
            
            if intensity < baseline:
                continue
            
            if rt >= min_x and rt <= max_x:
                refined_x.append(rt)
                refined_y.append(intensity)
        
        log_SSE = -1
        if min_sse > 0:
            log_SSE = np.log10(min_sse)
            log_SSE = round(log_SSE, 5)
            #print(min_sse, " log10:", log_SSE)
            
        #print("Refined Min:max=>", refined_x, ", ", refined_y)
        
        return refined_x, refined_y, log_SSE
    ######################## 2Demension ########################
    def doBoltzmann_refine(self, data_x, data_y, center, org_sigma, height, baseline):
        #print("ORG Sigma:", org_sigma)
        samplesize = 10000
        max_intensity = height
        loc = 0 #my guess
        
        sigma_arr = np.linspace(0, 0.01, 10)
        
        optimized_sigma= -1
        sse_list = list()
        min_sse = -1
        min_sse_rt = -1
        is_first= True
        for sigma_idx in range (0, len(sigma_arr)):
            sigma = org_sigma - sigma_arr[sigma_idx]
            a_value = np.sqrt((sigma**2 * math.pi)/(3*math.pi - 8)) #calculated based on wiki description
            v_2d = maxwell.rvs(loc, a_value, size=samplesize) #array corresponding to 2D proper motion obtained from Hubbs
            mean, var, skew, kurt = maxwell.stats(moments='mvsk')
            maxx = np.linspace(min(v_2d), max(v_2d), samplesize)
            params = maxwell.fit(v_2d, floc=0)
            factor_y = max_intensity / np.max(maxwell.pdf(maxx, *params))
            new_y = maxwell.pdf(maxx, *params) * factor_y
            max_var = max(new_y)
            max_idx = -1
            for i in range(0, len(new_y)):
                if new_y[i] == max_var:
                    max_idx = i
                    break
            max_var_x = maxx[max_idx]
            factor_x = center - max_var_x
            new_x = maxx + factor_x
            
            # SSE 계산 
            shift_rt_arr = np.linspace(-0.1, 0.1, 10)
            for idx in range(0, len(shift_rt_arr)):
                sse_list.clear()
                shift_rt = shift_rt_arr[idx]
                
                shifted_rt_arr = self.rt_shift(data_x, shift_rt)
                
                sse_list = list()
                
                for i in range(0, len(shifted_rt_arr)):
                    x = shifted_rt_arr[i] - factor_x
                    
                    x_min = min(new_x)
                    x_max = max(new_x)
                    if shifted_rt_arr[i] < x_min or x > x_max:
                        continue
                    
                    _var = maxwell.pdf(x, *params)
                    y = data_y[i]
                    _var_weight = _var * factor_y
                    sse = (y - _var_weight) **2
                    sse_list.append(sse)

                if is_first:
                    min_sse = sum(sse_list)
                    min_sse_rt = shift_rt
                    optimized_sigma = sigma
                    is_first = False
                else:
                    if min_sse > sum(sse_list):
                        min_sse = sum(sse_list)
                        min_sse_rt = shift_rt
                        optimized_sigma = sigma
                #print("Sigma:", sigma, " SSE:", sum(sse_list))
        #print("Optimized Sigma:", optimized_sigma, " Minimum SSE : ", min_sse_rt)
        optimized_x = self.rt_shift(new_x, -min_sse_rt)
        a_value = np.sqrt((optimized_sigma**2 * math.pi)/(3*math.pi - 8)) #calculated based on wiki description
        v_2d = maxwell.rvs(loc, a_value, size=samplesize) #array corresponding to 2D proper motion obtained from Hubbs
        mean, var, skew, kurt = maxwell.stats(moments='mvsk')
        maxx = np.linspace(min(v_2d), max(v_2d), samplesize)
        params = maxwell.fit(v_2d, floc=0)
        factor_y = max_intensity / np.max(maxwell.pdf(maxx, *params))
        new_y = maxwell.pdf(maxx, *params) * factor_y
        max_var = max(new_y)
        max_idx = -1
        for i in range(0, len(new_y)):
            if new_y[i] == max_var:
                max_idx = i
                break
        max_var_x = maxx[max_idx]
        factor_x = center - max_var_x
        new_x = maxx + factor_x
        
        
        # 범위 내의 +- 0.3 구간출력 
        refined_x = list();
        refined_y = list();
        
        min_x = np.min(data_x) - 0.3
        max_x = np.max(data_x) + 0.3
        
        for i in range(0, len(new_x)):
            rt = new_x[i]
            intensity = new_y[i]
            
            if intensity < baseline:
                continue
            
            #if rt >= final_min_rt and rt <= final_max_rt:
            if rt >= min_x and rt <= max_x:            
                refined_x.append(rt)
                refined_y.append(intensity)
        
        log_SSE = -1
        if sse > 0:
            log_SSE = np.log10(sse)
            log_SSE = round(log_SSE, 5)
        
        return refined_x, refined_y, log_SSE
    
    def rt_shift(self, data_x, shift_rt):
        rt_list = list()
        for i in range(0, len(data_x)):
            x = data_x[i]
            x_ = x + (shift_rt)
            rt_list.append(x_)
        
        return rt_list
        
    def get_model_SIL(self):
        return self.model_SIL
    def get_model_Endo(self):
        return self.model_Endo
    
    def get_area_SIL(self):
        return self.area_SIL
    def get_area_Endo(self):
        return self.area_Endo
    def get_area_ratio(self):
        return self.area_ratio    
    
class MRMProfile:
    def __init__(self, peptide_seq):
        self.peptide_seq = peptide_seq
        self.transitionList = list()
        self.modeled_transitionList = list()
        self.min_rt = 999999999.0
        self.max_rt = 0
        self.is_modeling = 0 #모델링이 끝났다면 1로 세팅 
        self.avg_ratio = -1
        self.img_path = None
        self.target_transition_str = ""        
        
    def set_rep_area_SIL(self, area):
        self.rep_area_SIL = area;
    def set_rep_area_Endo(self, area):
        self.rep_area_Endo = area
    
    def get_rep_area_SIL(self):
        return self.rep_area_SIL
    def get_rep_area_Endo(self):
        return self.rep_area_Endo
    
    def set_all_ratio(self, ratio):
        self.all_ratio = ratio
    def get_all_ratio(self):
        return self.all_ratio
        
    def set_target_transition_str(self, target_transition_str):
        self.target_transition_str= target_transition_str
    def get_target_transition_str(self):
        return self.target_transition_str
    
    def set_avg_ratio(self, avg_ratio):
        self.avg_ratio = avg_ratio
    def get_avg_ratio(self):
        return self.avg_ratio
    
    def set_avg_SIL_Area(self, avg_area):
        self.avg_SIL_Area = avg_area
    def get_avg_SIL_Area(self):
        return self.avg_SIL_Area
    
    def set_avg_SIL_Area_Str(self, avg_area_str):
        self.avg_SIL_Area_Str = avg_area_str
    def get_avg_SIL_Area_Str(self):
        return self.avg_SIL_Area_Str
    
    def set_avg_Endo_Area_Str(self, avg_area_str):
        self.avg_Endo_Area_Str = avg_area_str
    def get_avg_Endo_Area_Str(self):
        return self.avg_Endo_Area_Str
    
    def set_avg_rt(self, avg_rt):
        self.avg_rt = avg_rt
    def get_avg_rt(self):
        return self.avg_rt
    
    def set_avg_Endo_Area(self, avg_area):
        self.avg_Endo_Area = avg_area
    def get_avg_Endo_Area(self):
        return self.avg_Endo_Area
    
    def set_is_modeling(self, is_modeling):
        self.is_modeling = is_modeling
    def get_is_modeling(self):
        return self.is_modeling
        
    def set_min_rt(self, min_rt):
        self.min_rt = min_rt
    def set_max_rt(self, max_rt):
        self.max_rt = max_rt
    
    def setPeptideSeq(self, peptide_seq):
        self.peptide_seq = peptide_seq
    def getPeptideSeq(self):
        return self.peptide_seq
    
    def setPrecursorMz(self, precursor_mz):
        self.precursor_mz = precursor_mz
    def getPrecursorMz(self):
        return self.precursor_mz
    
    def setRt(self, rt):
        self.rt = rt
    def getRt(self):
        return self.rt
    
    def setIntensitySIL(self, intensitySIL):
        self.intensitySIL = intensitySIL
    def getIntensitySIL(self):
        return self.intensitySIL
    
    def setIntensityEndo(self, intensityEndo):
        self.intensityEndo = intensityEndo
    def getIntensityEndo(self):
        return self.intensityEndo
    
    def setProfileSIL(self, profileSIL):
        self.profileSIL = profileSIL
    def setProfileEndo(self, profileEndo):
        self.profileEndo = profileEndo
    def getProfileSIL(self):
        return self.profileSIL
    def getProfileEndo(self):
        return self.profileEndo    
    
    def setRatio(self, ratio):
        self.ratio = ratio
    def getRatio(self):
        return self.ratio
    
    def setImgPath(self, path):
        self.img_path = path
    def getImgPath(self):
        return self.img_path
    
    def addTransition(self, transition):
        self.transitionList.append(transition)
    def getTransitionList(self):
        return self.transitionList
    def setTransitionList(self, transitionList):
        self.transitionList = transitionList
        
    def setModeled_TransitionList(self, transitionList):
        self.modeled_transitionList = transitionList;
    def getModeled_TransitionList(self):
        return self.modeled_transitionList;
    
    def getTransition(self, precursorMz, ionType):
        targetTransition = None
        for transition in self.transitionList:
            if transition.getPrecursorMz() == precursorMz and transition.getIonType() == ionType:
                targetTransition = transition
                break  
        
        if targetTransition == None:
            targetTransition = Transition(precursorMz, ionType)
            self.transitionList.append(targetTransition)
            return targetTransition
        else:
            return targetTransition
        
    def get_peak_width(self):
        targetTransition = None
        
        max_peak_width = 0.0
        
        for transition in self.transitionList:
            if transition != None:
                if len(transition.peak_width_list) > 0:
                    peak_width = max(transition.peak_width_list)
                    if max_peak_width < peak_width:
                        max_peak_width = peak_width
        
        return max_peak_width
        
    def getProfileSum(self):
        sumIntensitySIL = 0
        sumIntensityEndo = 0
        ratio = 0.0
        if len(self.transitionList) > 0:
            transition = self.transitionList[0]
            profileSIL = transition.getProfileSIL()
            profileEndo = transition.getProfileEndo()
            
            intensityListSIL = profileSIL.values()
            intensityListEndo = profileEndo.values()
            
            if len(intensityListSIL) > 0:
                sumIntensitySIL = sum(intensityListSIL)
                sumIntensityEndo = sum(intensityListEndo)
                
                '''
                for i in range(0, len(intensityListSIL)):
                    intensitySIL = intensityListSIL[i]
                    intensityEndo = intensityListEndo[i]
                    
                    sumIntensitySIL += intensitySIL
                    sumIntensityEndo += intensityEndo
                '''
                if sumIntensitySIL > 0:
                    ratio = sumIntensityEndo / sumIntensitySIL
        return sumIntensitySIL, sumIntensityEndo, ratio
            
    def do_modeling(self):
        sorted_transition_list = list()
        for tr in self.transitionList:
            gen_model = tr.fit_model()
            if gen_model == 0:
                continue
            
            model_SIL = tr.get_model_SIL()[tr.get_final_model_index_SIL()]
            min_model_x, max_model_x = model_SIL.get_model_range()
            model_Endo = tr.get_model_Endo()[tr.get_final_model_index_Endo()]
            #ratio = round(model_Endo.get_model_area() / model_SIL.get_model_area(),4)
            sil_model_area = model_SIL.get_model_area()
            endo_model_area = model_Endo.get_model_area_with_range(min_model_x, max_model_x)
            ratio = round(endo_model_area / model_SIL.get_model_area(),4)
            
            self.ratio_str += tr.getIonType() + ":" + str(ratio) + "\n"
            
            rt_list_SIL = model_SIL.get_x()
            if self.min_rt > min(rt_list_SIL):
                self.min_rt = min(rt_list_SIL)
            if self.max_rt < max(rt_list_SIL):
                self.max_rt = max(rt_list_SIL)
            
            #self.do_sort(sorted_transition_list, tr)
            sorted_transition_list.append(tr)
        self.transitionList = sorted_transition_list        
                
    def sort_by_area_custom(self, model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, m_factor):
        #print("No. transitions:", len(self.transitionList))
        sorted_transition_list = list()
        for tr in self.transitionList:
            
            profileSIL = tr.getProfileSIL()
            profileEndo = tr.getProfileEndo()
            
            rtListORG = list(profileSIL.keys())
            intensityListSIL = list(profileSIL.values())
            rtList, intensityListSIL = tr.check_range(rtListORG, intensityListSIL, rt_from, rt_to)
        
            intensityListEndo = list(profileEndo.values())        
            rtList, intensityListEndo = tr.check_range(rtListORG, intensityListEndo, rt_from, rt_to)
        
            m_rt_intensity_SIL = tr.getMovingAverage(rtList, intensityListSIL, m_factor)        
            m_rt_intensity_Endo = tr.getMovingAverage(rtList, intensityListEndo, m_factor)
        
            tr.setProfileSIL(m_rt_intensity_SIL)
            tr.setProfileEndo(m_rt_intensity_Endo)
            
            #print("fit modeling:", tr.getPrecursorMz(), "->",tr.getTransitionMz())
            tr.fit_model_custom(model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, m_factor)
            
            sorted_transition_list.append(tr)
        
        self.transitionList = sorted_transition_list
        
        for tr in self.transitionList:
            
            model_SIL = tr.get_model_SIL()[tr.get_final_model_index_SIL()]
            model_Endo = tr.get_model_Endo()[tr.get_final_model_index_Endo()]
            ratio = round(model_Endo.get_model_area() / model_SIL.get_model_area(),4)
            
            self.ratio_str += tr.getIonType() + ":" + str(ratio) + "\n"
            
            rt_list_SIL = model_SIL.get_x()
            if self.min_rt > min(rt_list_SIL):
                self.min_rt = min(rt_list_SIL)
            if self.max_rt < max(rt_list_SIL):
                self.max_rt = max(rt_list_SIL)                
            
    def do_sort(self, sorted_list, transition):
        index = -1
        for i in range(0, len(sorted_list)):
            tr = sorted_list[i]
            if tr.get_area_SIL() < transition.get_area_SIL():
                index = i
                break
        
        if index == -1:
            sorted_list.append(transition)
        else:
            sorted_list.insert(index, transition)
            
    def get_ratio_str(self):
        return self.ratio_str
    
    def get_min_rt(self):
        return float(self.min_rt) - 0.15
    def get_max_rt(self):
        return float(self.max_rt) + 0.15   

class CanvasWidget(QWidget):

    def __init__(self, parent, peptide_seq, transition, transition_index):        
        super(CanvasWidget, self).__init__(parent)
        self.parent = parent
        self.peptide_seq = peptide_seq
        self.transition = transition
        self.transition_type = transition.getIonType()
        self.transition_index = transition_index
        self.profileSIL = transition.getProfileSIL()
        self.profileEndo = transition.getProfileEndo()
        self.rtList = list(self.profileSIL.keys())
        self.intensityListSIL = list(self.profileSIL.values())
        self.intensityListEndo = list(self.profileEndo.values())
        #self.min_x = transition.get_rt_min_for_view() - 0.3
        #self.max_x = transition.get_rt_max_for_view() + 0.3
        self.min_x = transition.get_rt_min_for_view()
        self.max_x = transition.get_rt_max_for_view()
        self.sse_SIL = -1
        self.sse_Endo = -1
        self.model_area_SIL = 0
        self.model_area_Endo = 0
        self.ratio = 0
        self.rt_min_for_view = self.parent.rt_min_for_view
        self.rt_max_for_view = self.parent.rt_max_for_view
        
        self.color_list = [ "red", "green", "orange", "purple", "brown", "pink", "gray", "#808000", "#00FFFF", "#8B0000", "#006400", "#FF8C00", "#EE82EE", "#8B4513", "#FF00FF", "#008B8B", "#8B008B" ]
        self.initUI()
        
    def isHigherSIL(self):
        maxSIL = max(self.intensityListSIL)
        maxEndo = max(self.intensityListEndo)
        
        if maxSIL > maxEndo:
            return True
        else:
            return False
        
    def get_transition(self):
        return self.transition
        
    def initUI(self):
            
        _layout = QVBoxLayout()
        _groupBox = QGroupBox()
        _groupBox.setFixedSize(650, 370)
        
        _boxlayout = QVBoxLayout()
        _1_box_layout = QHBoxLayout()
        _2_box_layout = QHBoxLayout()
        _3_box_layout = QHBoxLayout()
        _4_box_layout = QHBoxLayout()
        _blank_layout = QHBoxLayout()
        
        self.fig = plt.Figure(figsize=(4,5), dpi=(80))
        plt.subplots_adjust()
        self._canvas = FigureCanvas(self.fig)
        self._canvas.setFixedSize(450, 250)
        self.ax = self._canvas.figure.subplots()
        
        ############### Graph 출력 #################
        self.title = self.transition.getIonType() + "\n(" + "{:.3f}".format(float(self.transition.getPrecursorMz())) + "Th -> " + "{:.3f}".format(float(self.transition.getTransitionMz())) + "Th)\n"
        #self.ax.set_title(self.title)
        
        self.ax.set_xlim(self.min_x, self.max_x)
        self.new_rt_list_SIL, self.new_intensity_list_SIL = self.check_range(self.rtList, self.intensityListSIL, self.min_x, self.max_x)
        self.new_rt_list_Endo, self.new_intensity_list_Endo = self.check_range(self.rtList, self.intensityListEndo, self.min_x, self.max_x)
        self.ax.scatter(self.new_rt_list_SIL, self.new_intensity_list_SIL, marker=".", s=10, color=self.color_list[self.transition_index])
        self.ax.scatter(self.new_rt_list_SIL, self.new_intensity_list_Endo, marker=".", s=10, color="blue")
        
        
        self.model_SIL = None
        self.model_SIL_cnt = 0

        self.model_Endo = None
        self.model_Endo_cnt = 0
        if len(self.transition.get_model_Endo()) > 0:
            model_combo_index = self.transition.get_final_model_index_Endo()
            if model_combo_index >= len(self.transition.get_model_Endo()):
                model_combo_index = 0
            self.model_Endo = self.transition.get_model_Endo()[model_combo_index]
            self.model_Endo_cnt = len(self.transition.get_model_Endo())
            self.new_x_Endo = self.model_Endo.get_x()
            self.new_y_Endo = self.model_Endo.get_y()
            self.ax.plot(self.new_x_Endo, self.new_y_Endo, color="blue")
            self.sse_Endo = self.model_Endo.get_sse()
            self.model_area_Endo = self.model_Endo.get_model_area()        
            
        
        #print("Model SIL:", self.transition.get_model_SIL(), ":", self.transition.get_model_Endo())
        if len(self.transition.get_model_SIL()) > 0:
            model_combo_index = self.transition.get_final_model_index_SIL()
            if model_combo_index >= len(self.transition.get_model_SIL()):
                model_combo_index = 0
            self.model_SIL = self.transition.get_model_SIL()[model_combo_index]
            self.model_SIL_cnt = len(self.transition.get_model_SIL())
            self.new_x_SIL = self.model_SIL.get_x()
            self.new_y_SIL = self.model_SIL.get_y()
            self.ax.plot(self.new_x_SIL, self.new_y_SIL, color=self.color_list[self.transition_index])
            self.sse_SIL = self.model_SIL.get_sse()
            self.model_area_SIL = self.model_SIL.get_model_area()
            self.ratio = round(self.model_area_Endo / self.model_area_SIL,4)
            
        ####### 서브 plot ######
        try:
            if self.isHigherSIL():
                axins = inset_axes(self.ax, width='30%', height='30%', loc='upper right', borderpad=1.0)
                axins.scatter(self.new_rt_list_Endo, self.new_intensity_list_Endo, marker=".", s=10, color="blue")
                if max(self.new_intensity_list_Endo) > 0:
                    axins.plot(self.new_x_Endo, self.new_y_Endo, label='Graph 2', color='blue')
                    #axins.set_xlabel('')
                    #axins.set_ylabel('', color='r')
                axins.tick_params(axis='y', labelcolor='r')
            else:
                axins = inset_axes(self.ax, width='30%', height='30%', loc='upper right', borderpad=1.0)
                axins.scatter(self.new_rt_list_SIL, self.new_intensity_list_SIL, marker=".", s=10, color=self.color_list[self.transition_index])
                if max(self.new_y_SIL) > 0:
                    axins.plot(self.new_x_SIL, self.new_y_SIL, label='Graph 2', color=self.color_list[self.transition_index])
                #axins.set_xlabel('')
                #axins.set_ylabel('', color='r')
                axins.tick_params(axis='y', labelcolor='r')
        except AttributeError as ae:   
            pass
        
        
        #x_label = "ratio:" + str(ratio) + "\nSSE(SIL):" + str(sse_SIL) + "\nSSE(Endo):" + str(sse_Endo)
        self.x_label = "SSE(SIL):" + str(self.sse_SIL) + "\nSSE(Endo):" + str(self.sse_Endo)
        #self.ax.set_xlabel(self.x_label)
        
        
        self._canvas.draw()
        _transition_label = QLabel(self.transition_type)
        #_label = QLabel("Replicate #1")
        #_label.setStyleSheet("QLabel { color:#ff0000;}")
        
        self._cb = QCheckBox("Use Transition")
        _model_label = QLabel("No. Models   SIL")            
        self._combo_SIL = QComboBox()
        self._combo_SIL.addItem("1")
        self._combo_SIL.addItem("2")
        self._combo_SIL.setCurrentIndex(self.model_SIL_cnt - 1)
        _model_SIL_idx_label = QLabel("Model idx     SIL")
        self._combo_SIL_idx = QComboBox()
        for idx in range(0, self.model_SIL_cnt):
            self._combo_SIL_idx.addItem(str(idx))            
        self._combo_SIL_idx.setCurrentIndex(self.transition.get_final_model_index_SIL())
        
        _model_label_Endo = QLabel("         Endo")
        self._combo_Endo = QComboBox()
        self._combo_Endo.addItem("1")
        self._combo_Endo.addItem("2")
        self._combo_Endo.setCurrentIndex(self.model_Endo_cnt - 1)
        
        _model_Endo_idx_label = QLabel("         Endo")
        self._combo_Endo_idx = QComboBox()
        for idx in range(0, self.model_Endo_cnt):
            self._combo_Endo_idx.addItem(str(idx))
        self._combo_Endo_idx.setCurrentIndex(self.transition.get_final_model_index_Endo())
        
        self._apply_button = QPushButton("Apply")
        self._refit_button = QPushButton("Fit")
        
        _blanks_label = QLabel("     ")
        _blanks_label2 = QLabel("  ")
        
        #### canvas 옆 테이블 ####
        self.transition_info_table = QTableWidget()
        self.transition_info_table.verticalScrollBar().setVisible(False)
        self.transition_info_table.setFixedSize(180, 210)
        self.transition_info_table.horizontalHeader().setVisible(False)
        self.transition_info_table.setRowCount(7)
        self.transition_info_table.setColumnCount(1)
        column_headers = ['IonType', 'm/z', 'SSE(SIL)', 'SSE(Endo)  ', 'Ratio', 'Replicate #1', 'Replicate #2']
        self.transition_info_table.setVerticalHeaderLabels(column_headers)
        self.transition_info_table.verticalHeader().setStyleSheet("QHeaderView::section {background-color:#404040;color:#FFFFFF;}")
        #self.transition_info_table.insertColumn(1)
        #mz_str = str(self.transition.getPrecursorMz()) + "Th ->\n " + str(self.transition.getTransitionMz())
        mz_str = str(round(float(self.transition.getTransitionMz()),4)) + "Th"
        self.transition_info_table.setItem(0, 0, QTableWidgetItem(self.transition.getIonType()))
        self.transition_info_table.setItem(1, 0, QTableWidgetItem(mz_str))
        self.transition_info_table.setItem(2, 0, QTableWidgetItem(str(self.sse_SIL)))
        self.transition_info_table.setItem(3, 0, QTableWidgetItem(str(self.sse_Endo)))
        self.transition_info_table.setItem(4, 0, QTableWidgetItem(str(self.ratio)))
        font = self.transition_info_table.font()
        font.setPointSize(10)
        self.transition_info_table.setFont(font)
        
        is_rep1 = 0
        if self.peptide_seq in self.parent.replicate1_check_table: 
            is_rep1 = self.parent.check_replicate1(self.peptide_seq, self.transition.getIonType())
        
        if is_rep1 == 0:
            self.transition_info_table.setItem(5, 0, QTableWidgetItem("No"))
        else:
            self.transition_info_table.setItem(5, 0, QTableWidgetItem("Yes"))
            
        is_rep2 = 0
        if self.peptide_seq in self.parent.replicate2_check_table: 
            is_rep2 = self.parent.check_replicate2(self.peptide_seq, self.transition.getIonType())
        if is_rep2 == 0:
            self.transition_info_table.setItem(6, 0, QTableWidgetItem("No"))
        else:
            self.transition_info_table.setItem(6, 0, QTableWidgetItem("Yes"))
        
        #_boxlayout.addWidget(self._canvas)
        _1_box_layout.addWidget(self._canvas)
        _1_box_layout.addWidget(self.transition_info_table)
        
        _2_box_layout.addWidget(self._cb)

        _3_box_layout.addWidget(_model_label)
        _3_box_layout.addWidget(self._combo_SIL)
        _3_box_layout.addWidget(_model_label_Endo)
        _3_box_layout.addWidget(self._combo_Endo)
        _3_box_layout.addWidget(_blanks_label2)
        _3_box_layout.addWidget(self._refit_button)
        
        _4_box_layout.addWidget(_model_SIL_idx_label)
        _4_box_layout.addWidget(self._combo_SIL_idx)
        _4_box_layout.addWidget(_model_Endo_idx_label)
        _4_box_layout.addWidget(self._combo_Endo_idx)
        _4_box_layout.addWidget(_blanks_label2)
        _4_box_layout.addWidget(self._apply_button)
        
        _blank_layout.addWidget(_blanks_label)
        
        _boxlayout.addLayout(_1_box_layout)
        _boxlayout.addLayout(_blank_layout)
        _boxlayout.addLayout(_2_box_layout)
        _boxlayout.addLayout(_3_box_layout)
        _boxlayout.addLayout(_4_box_layout)
        
        _groupBox.setLayout(_boxlayout)
        _layout.addWidget(_groupBox)
        
        
        if self.transition.get_use_transition() == 1:
            self._cb.toggle()
            self.parent.add_target_transition(self.transition)
            
        self._refit_button.clicked.connect(self.refit)
        self._apply_button.clicked.connect(self.change_model)
        self._canvas.mpl_connect('button_press_event', self.canvas_clicked)
        self._cb.stateChanged.connect(self.useTransition)
        self.setLayout(_layout)
    
    def refit(self):
        
        if self._combo_SIL.currentText() == "":
            self._combo_SIL.setCurrentIndex(0)
            
        if self._combo_Endo.currentText() == "":
            self._combo_Endo.setCurrentIndex(0)
            
        model_cnt_SIL = int(self._combo_SIL.currentText())
        model_cnt_Endo = int(self._combo_Endo.currentText())
        self.parent.refit(self.transition, model_cnt_SIL, model_cnt_Endo)
        
    def change_model(self):
        
        if self._combo_SIL_idx.currentText() == "":
            self._combo_SIL_idx.setCurrentIndex(0)
        
        if self._combo_Endo_idx.currentText() == "":
            self._combo_Endo_idx.setCurrentIndex(0)
        
        model_idx_SIL = int(self._combo_SIL_idx.currentText())
        model_idx_Endo = int(self._combo_Endo_idx.currentText())
        
        try:
            self.parent.change_model(self.transition, model_idx_SIL, model_idx_Endo)
        except TypeError as ae:   
                self.show_messageBox("오류", "예기치 못한 오류가 발생 했습니다\n다시 시도 하세요")
        
    def canvas_clicked(self, event):
        self._cb.toggle()
        
    def useTransition(self):
        #print("call useTransition() ...")
        if self._cb.isChecked() == True:
            self.parent.add_target_transition(self.transition)
        else:
            self.parent.remove_target_transition(self.transition)
            
    def check_range(self, rtList, intensityList, min_x, max_x):
        
        #print("check_range:", min_x, "~", max_x)
        
        new_rt_list = list()
        new_intensity_list = list()
        
        for i in range(0, len(rtList)):
            rt = rtList[i]
            if rt < min_x or rt > max_x:
                continue
            new_rt_list.append(rt)
            new_intensity_list.append(intensityList[i])
        
        return new_rt_list, new_intensity_list

class MainWindow(QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.initUI()
        self.target_transition_list = list()
        #self.threadclass = ThreadClass()
        self.color_list = [ "red", "green", "orange", "purple", "brown", "pink", "gray", "#808000", "#00FFFF", "#8B0000", "#006400", "#FF8C00", "#EE82EE", "#8B4513", "#FF00FF", "#008B8B", "#8B008B" ]
        self.replicate1_check_table = {}
        self.replicate2_check_table = {}
        self.saveCnt = 0
        self.step = 100
        self.rt_min_for_view = 0
        self.rt_max_for_view = 0
        
    def initUI(self):
        
        #self.setGeometry(50, 70, 1600,1200)
        self.setWindowTitle("Wideband PRM")
        self.statusbar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusbar)
        self.progressbar = QProgressBar()
        self.statusbar.addPermanentWidget(self.progressbar)
        
        mainWidget = QWidget()
        main_layout = QHBoxLayout()
        mainWidget.setLayout(main_layout)
        
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        
        #################### replicate load button ######################################
        _replicate_groupbox = QGroupBox()
        _replicate_groupbox_layout = QVBoxLayout()
        _rep_1_layout = QHBoxLayout()
        _rep_2_layout = QVBoxLayout()
        _replicate_groupbox_layout.addLayout(_rep_1_layout)
        _replicate_groupbox_layout.addLayout(_rep_2_layout)
        
        _replicate_label = QLabel("Replicate Results")
        self._rep1_button = QPushButton("Replicate #1")
        self._rep2_button = QPushButton("Replicate #2")
        
        _rep_1_layout.addWidget(_replicate_label)
        _rep_1_layout.addWidget(self._rep1_button)
        _rep_1_layout.addWidget(self._rep2_button)
        
        self._selected_rep1_label = QLabel("No Replicate #1")
        self._selected_rep2_label = QLabel("No Replicate #2")
        _rep_2_layout.addWidget(self._selected_rep1_label)
        _rep_2_layout.addWidget(self._selected_rep2_label)
        
        _replicate_groupbox.setLayout(_replicate_groupbox_layout)
        left_layout.addWidget(_replicate_groupbox)
        #######################################################################
        
        #################### load button ######################################
        _load_groupbox = QGroupBox()
        _load_groupbox_layout = QHBoxLayout()
        
        _load_label = QLabel("File")
        self._load_tf = QLineEdit()
        self._load_tf.setFixedWidth(580)
        self._load_button = QPushButton("Load")
        _load_groupbox_layout.addWidget(_load_label)
        _load_groupbox_layout.addWidget(self._load_tf)
        _load_groupbox_layout.addWidget(self._load_button)
        _load_groupbox.setLayout(_load_groupbox_layout)
        left_layout.addWidget(_load_groupbox)
        #######################################################################
        
        #################### Parameter ######################################
        _param_groupbox = QGroupBox()
        _param_groupbox_layout = QHBoxLayout()
        _param_label = QLabel("RT")
        self._rt_from_tf = QLineEdit()
        self._rt_from_tf.setFixedWidth(80)
        _till_label = QLabel("    ~")
        self._rt_to_tf = QLineEdit()
        self._rt_to_tf.setFixedWidth(80)
        
        _global_SIL_label = QLabel("        No. Model (SIL)")
        self._global_SIL_combo = QComboBox()
        self._global_SIL_combo.addItem("1")
        self._global_SIL_combo.addItem("2")
        
        _global_Endo_label = QLabel("   No. Model (Endo)")
        self._global_Endo_combo = QComboBox()
        self._global_Endo_combo.addItem("1")
        self._global_Endo_combo.addItem("2")
        
        _blank_label = QLabel("             ")
        
        self._fit_button = QPushButton("Fit All")
        _param_groupbox_layout.addWidget(_param_label)
        _param_groupbox_layout.addWidget(self._rt_from_tf)
        _param_groupbox_layout.addWidget(_till_label)
        _param_groupbox_layout.addWidget(self._rt_to_tf)
        _param_groupbox_layout.addWidget(_global_SIL_label)
        _param_groupbox_layout.addWidget(self._global_SIL_combo)
        _param_groupbox_layout.addWidget(_global_Endo_label)
        _param_groupbox_layout.addWidget(self._global_Endo_combo)
        _param_groupbox_layout.addWidget(_blank_label)
        _param_groupbox_layout.addWidget(self._fit_button)
        _param_groupbox.setLayout(_param_groupbox_layout)
        left_layout.addWidget(_param_groupbox)
        #######################################################################
        
        canvas_layout = QHBoxLayout()
        _1_canvas_layout = QVBoxLayout()
        _2_canvas_layout = QVBoxLayout()
        canvas_layout.addLayout(_1_canvas_layout)
        canvas_layout.addLayout(_2_canvas_layout)
        
        ########## Canvas 부분 구성 #############
        widget = QWidget()
        layout = QVBoxLayout()
        widget.setLayout(layout)
        widget.setFixedSize(730, 810)
        _1_canvas_layout.addWidget(widget)
        
        groupBox = QtWidgets.QGroupBox(widget)
        groupBox.setGeometry(QtCore.QRect(0, 0, 780, 800))
        groupBox.setObjectName("groupBox")
        
        self.scroll = QScrollArea(groupBox)
        self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scroll.verticalScrollBar().setEnabled(True)
        self.scroll.horizontalScrollBar().setEnabled(False)
        self.scroll.setWidgetResizable(True)
        self.scroll.setFixedSize(730, 800)
        
        canvasWidget = QWidget(self.scroll)
        self._layout = QVBoxLayout(canvasWidget)
        ############### 요기에 각 transition별로 그래프가 출력되면 된다 ##############
        

        self.scroll.setWidget(canvasWidget)
        
        ########## Canvas 우측 부분 버튼 구성 #####
        self.save_ratio_button = QPushButton("\n\n\n\n\n\n\n\n\n\n\n\n\n\nSave Ratio\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
        #self.save_ratio_button = QPushButton("Save Ratio")
        self.save_ratio_button.setFont(QtGui.QFont("Arial", 13, QtGui.QFont.Bold))
        
        self.save_ratio_button.setStyleSheet("color:white;background-color: #1E90FF;border-width: 3px;")
        
        _2_canvas_layout.addWidget(self.save_ratio_button)        
        
        #left_layout.addWidget(widget)
        left_layout.addLayout(canvas_layout)
        
        #################### save button ######################################
        _save_groupbox = QGroupBox()
        #_save_groupbox.setGeometry(QtCore.QRect(10, 600, 645, 768))
        #_save_groupbox.setFixedSize(645, 50)
        _save_groupbox_layout = QHBoxLayout()
        
        _save_label = QLabel("Result ")
        self._save_tf = QLineEdit()
        _fileChoose_button = QPushButton("...")
        self._save_button = QPushButton("Save")
        _save_groupbox_layout.addWidget(_save_label)
        _save_groupbox_layout.addWidget(self._save_tf)
        #_save_groupbox_layout.addWidget(_fileChoose_button)
        _save_groupbox_layout.addWidget(self._save_button)
        _save_groupbox.setLayout(_save_groupbox_layout)
        left_layout.addWidget(_save_groupbox)
        #######################################################################
        
        ######################### peptide list table ##########################
        self.peptide_list_table = QTableWidget()
        self.peptide_list_table.setColumnCount(3)
        column_headers = ['               Peptide               ', ' Precursor m/z ', '      Ratio      ']
        self.peptide_list_table.setHorizontalHeaderLabels(column_headers)
        self.peptide_list_table.horizontalHeader().setStyleSheet("QHeaderView::section {background-color:#404040;color:#FFFFFF;}")
        right_layout.addWidget(self.peptide_list_table)
        #######################################################################
        
        ##################### 이벤트 등록 ############################
        self._rep1_button.clicked.connect(self.load_replicate1)
        self._rep2_button.clicked.connect(self.load_replicate2)
        self._fit_button.clicked.connect(self.fit_all)
        self.save_ratio_button.clicked.connect(self.save_ratio)
        self.peptide_list_table.itemClicked.connect(self.showCurrentRow)        
        self._load_button.clicked.connect(self.load_input)
        self._save_button.clicked.connect(self.save)
        
        self.setFixedSize(1350,1250)
        self.setCentralWidget(mainWidget)
        self.currentWidgets = list()
        self.show()
        
        ############## 자동화 #############
        #self._load_tf.setText("M:/Python/Test/WBPRM_2X_PDAC030_3ug_Cycl1sec_RT6min_WatchRes15_120min_rep1_C1_021622.raw")
        #self.input_file_name = self._load_tf.text().strip()
        #self.load_file()
        
    
    def add_target_transition(self, transition):
        if transition not in self.target_transition_list:
            self.target_transition_list.append(transition)        
    
    def remove_target_transition(self, transition_type):
        if transition_type in self.target_transition_list:
            self.target_transition_list.remove(transition_type)
    
    def show_msg(self, msg):
        self.statusbar.showMessage(msg)

    def load_replicate1(self):
        self.rep1_file_name = QtWidgets.QFileDialog.getOpenFileName(self, "file select Replicate1")[0]
        if self.rep1_file_name == "":
            return
        self._selected_rep1_label.setText("Replicate1 : " + self.rep1_file_name)
        df = pd.read_excel(self.rep1_file_name)
        rep1_peptide_list = df.iloc[:, 0]
        rep1_transition_list = df.iloc[:, 2]
        
        for i in range(0, len(rep1_peptide_list)):
            peptide = rep1_peptide_list[i]
            transitions = rep1_transition_list[i].split("\n")
            self.replicate1_check_table[peptide] = transitions
            #print(peptide, ":", transitions)
            
    def load_replicate2(self):
        self.rep2_file_name = QtWidgets.QFileDialog.getOpenFileName(self, "file select Replicate2")[0]
        if self.rep2_file_name == "":
            return
        self._selected_rep2_label.setText("Replicate2 : " + self.rep2_file_name)
        df = pd.read_excel(self.rep2_file_name)
        rep2_peptide_list = df.iloc[:, 0]
        rep2_transition_list = df.iloc[:, 2]
        
        for i in range(0, len(rep2_peptide_list)):
            peptide = rep2_peptide_list[i]
            transitions = rep2_transition_list[i].split("\n")
            self.replicate2_check_table[peptide] = transitions
            #print(peptide, ":", transitions)
            
    def check_replicate1(self, peptide, this_transition):
        rep1_transitions = self.replicate1_check_table[peptide]
        is_rep1 = 0
        for transition in rep1_transitions:
            if transition == this_transition:
                is_rep1 = 1
                break
        return is_rep1
    
    def check_replicate2(self, peptide, this_transition):
        rep2_transitions = self.replicate2_check_table[peptide]
        is_rep2 = 0
        for transition in rep2_transitions:
            if transition == this_transition:
                is_rep2 = 1
                break
        return is_rep2   
    
    
    def load_input(self):
        self.file_name = QtWidgets.QFileDialog.getOpenFileName(self, "file select", "./", "mzML File(*.mzML)")[0]
        if self.file_name == "":
            return
        #print(self.file_name)
        self._load_tf.setText(self.file_name)
        self.input_file_name = self._load_tf.text().strip()
        self.load_file()
        
    def save(self):
        self.output_file_name = QtWidgets.QFileDialog.getSaveFileName(self, "file select", "./","Excel(*.xlsx)")[0]
        #print("output:", self.output_file_name)
        if self.output_file_name == "":
            return
        
        self._save_tf.setText(self.output_file_name)
        try:
            self.make_excel(self.output_file_name)
        except xlsxwriter.exceptions.FileCreateError as ae:  
            self.show_messageBox("에러", "결과파일을 생성할 수 없습니다. 사용중인지 확인하세요")
            return
        
    def load_file(self):
        
        args = [
            "-TARGET1", "target1.csv",
            "-TARGET2", "target2.csv",
            "-TARGET3", "target3.csv",
            "-MZML", self.input_file_name,
            "-PEAKCOUNT", "100",
            "-TOLERANCE", "10.0"
        ]
        
        self.profileList = {}
        
        if self.input_file_name.endswith(".mzML"):
            parser = MzMLParser(self.input_file_name)
            spectrum_list, spectrum_table, scan_rt_table = parser.parse()
            print("No. spectrra : ", len(spectrum_list))
            
            wideband_prm = WidebandPRM(args, spectrum_list, spectrum_table, scan_rt_table)
            self.peak_tracing_file = wideband_prm.execute()
            
            self.loadPRMTarget()
            
            print("Peack Tracing : ", self.peak_tracing_file)
        
        # 기존의 self.peptide_list_table의 row가 있었다면 다 제거한다 
        rowCnt = self.peptide_list_table.rowCount()
        #print(rowCnt, " rows")
        if rowCnt > 0:
            for i  in range (0, rowCnt):
                self.peptide_list_table.removeRow(0)
        
        peptide_seq = ""
        precursor_mz = 0.0
        
        f = open(self.peak_tracing_file, 'r')
        line = f.readline()  # first comment line
        
        previous_peptide_seq = ""
        
        mrmProfile = None
        
        target_peptide_seq = ""
        
        while True:
            line = f.readline()    
            if line == "":
                break
            
            cols = line.split('\t')
            peptide_seq = cols[0]
            precursor_mz = cols[1]
            scanNum = cols[2]
            ion_type = cols[3]
            transitionMz = cols[4]
            rt = float(cols[5])
            intensitySIL = float(cols[6])
            intensityEndo = float(cols[7])
            ratio = cols[8]
            
            if peptide_seq in self.profileList:
                mrmProfile = self.profileList[peptide_seq]
            else:
                mrmProfile = MRMProfile(peptide_seq)
                mrmProfile.setPrecursorMz(precursor_mz)
                self.profileList[peptide_seq] = mrmProfile
                self.targetList[peptide_seq] = "OK"
                
            if intensitySIL == 0 and intensityEndo == 0:
                continue                
              
            transition = mrmProfile.getTransition(precursor_mz, ion_type)
            transition.setRt(rt)
            transition.setTransitionMz(transitionMz)
            transition.setScanNum(scanNum)
            transition.getProfileSIL()[rt] = intensitySIL
            transition.getProfileEndo()[rt] = intensityEndo
            previous_peptide_seq = peptide_seq
                    
        self.show_msg("No. Peptides : " + str(len(self.profileList)))
        self.view_peptide_list()
        
    def loadPRMTarget(self):
        
        dirURL = os.path.realpath(self.input_file_name)
        dirURL = os.path.abspath(os.path.join(dirURL, os.pardir))
        
        self.target1_URL = dirURL + "/Target1.csv"
        self.target2_URL = dirURL + "/Target2.csv"        
        
        # Target 파일 로딩
        self.targetList = {}
        f = open(self.target1_URL, 'r')
        line = f.readline()  # first comment line
        while True:
            line = f.readline()    
            if line == "":
                break
            
            cols = line.split(',')
            peptide_seq = cols[0]
            self.targetList[peptide_seq] = "NA"
        
        f = open(self.target2_URL, 'r')
        line = f.readline()  # first comment line
        while True:
            line = f.readline()    
            if line == "":
                break
            
            cols = line.split(',')
            peptide_seq = cols[0]
            self.targetList[peptide_seq] = "NA"

        print("No. targets : ", len(self.targetList))
        
    def view_peptide_list(self):
        #peptideList = self.profileList.keys()
        peptideList = self.targetList.keys()
        rowNum = 0
        for peptide_seq in peptideList:
            
            if peptide_seq in self.profileList:
                mrmProfile = self.profileList[peptide_seq]
                precursorMz = mrmProfile.getPrecursorMz()
                self.peptide_list_table.insertRow(rowNum)
                self.peptide_list_table.setItem(rowNum, 0, QTableWidgetItem(peptide_seq))
                self.peptide_list_table.setItem(rowNum, 1, QTableWidgetItem(precursorMz))
                rowNum = rowNum + 1
            
            else:
                precursorMz = "NA"
                self.peptide_list_table.insertRow(rowNum)
                self.peptide_list_table.setItem(rowNum, 0, QTableWidgetItem(peptide_seq))
                self.peptide_list_table.setItem(rowNum, 1, QTableWidgetItem(precursorMz))
                self.peptide_list_table.setItem(rowNum, 2, QTableWidgetItem(""))
                self.peptide_list_table.item(rowNum, 0).setBackground(QtGui.QColor(255,0,0))
                self.peptide_list_table.item(rowNum, 1).setBackground(QtGui.QColor(255,0,0))
                self.peptide_list_table.item(rowNum, 2).setBackground(QtGui.QColor(255,0,0))
                rowNum = rowNum + 1
        self.peptide_list_table.resizeColumnsToContents()
        self.peptide_list_table.resizeRowsToContents()
        
    ######################### Fit All 버튼 클릭시 ####################    
    def fit_all(self):
        
        self.rt_min_for_view = self._rt_from_tf.text()
        self.rt_max_for_view = self._rt_to_tf.text()
        
        print("call refit : ", self.rt_min_for_view, " ~ ", self.rt_max_for_view)
        
        self.progressbar.setValue(0)
        if len(self.currentWidgets) > 0:
            for w in self.currentWidgets:
                self._layout.removeWidget(w)
        
        #self.threadclass.start()
        self.show_msg(self.peptide_seq + " processing ...")
        
        if self._rt_from_tf.text().strip() == "" or self._rt_to_tf.text().strip() == "":
            self.show_messageBox("Info", "required RT condition!")
            return
        
        rt_from = float(self._rt_from_tf.text().strip())
        rt_to   = float(self._rt_to_tf.text().strip())
        model_cnt_SIL = int(self._global_SIL_combo.currentText())
        model_cnt_Endo = int(self._global_Endo_combo.currentText())       
        
        color_code = 0
        
        transitionList = self.mrmProfile.getTransitionList()
        new_transition_list = list()
        
        min_rt = 99999999
        max_rt = 0
        transitionIndex = 0
        for transition in transitionList:
            transition_type = transition.getIonType()
            transitionIndex += 1
            self.show_msg(transition_type + " 처리 중 ...")
            # progressbar 업데이트 
            self.step = int(transitionIndex / len(transitionList) * 100)
            self.progressbar.setValue(self.step)
            
            _val = transition.checkEndoPriority()
            if _val == 1:
                gen_model = transition.fit_model_custom_reverse(model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, 1)
                if gen_model == -1:
                    transitionIndex += 1
                    continue                
            else:
                gen_model = transition.fit_model_custom(model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, 1)
                if gen_model == -1:
                    transitionIndex += 1
                    continue
            
            model_SIL = transition.get_model_SIL()[0]
            model_Endo = transition.get_model_Endo()[0]
            ratio = round(model_Endo.get_model_area() / model_SIL.get_model_area(),4)
            rt_list_SIL = model_SIL.get_x()
            if max(rt_list_SIL) - min(rt_list_SIL) < 10:  # 2분 이상 인것들은 구간 측정에서 제외 
                if min_rt > min(rt_list_SIL):
                    min_rt = min(rt_list_SIL)
                if max_rt < max(rt_list_SIL):
                    max_rt = max(rt_list_SIL)            

            new_transition_list.append(transition)
        
        #print("No. new_transition_list:", len(new_transition_list))
        
        if min_rt != 99999999:
            self.mrmProfile.set_min_rt(self.rt_min_for_view)
            self.mrmProfile.set_max_rt(self.rt_max_for_view)
            self.mrmProfile.setModeled_TransitionList(new_transition_list)
            self._rt_from_tf.setText(str(round(min_rt,2)))
            self._rt_to_tf.setText(str(round(max_rt,2)))
            #self.mrmProfile.set_is_modeling(1)
        
        transitionIndex = 1
        transition_list = self.mrmProfile.getModeled_TransitionList()
        #print("No. modeled_transition_list:", len(transition_list))
        
        if len(transition_list) == 0:
            self.show_messageBox("정보", "생성된 모델이 없습니다. 조건을 변경하여 다시 시도하세요")
            return
        
        for transition in transition_list:
            transition.set_rt_min_for_view(self.mrmProfile.get_min_rt())
            transition.set_rt_max_for_view(self.mrmProfile.get_max_rt())
            
            # 그래프 구성 
            try:
                cWidget = CanvasWidget(self, self.mrmProfile.getPeptideSeq(), transition, transitionIndex)    
                self._layout.addWidget(cWidget)
                self.currentWidgets.append(cWidget)
                transitionIndex += 1             
            except AttributeError as ae:   
                print(ae)
        
        #self.statusbar.removeWidget(self.progressbar)
        
        self.show_msg(self.mrmProfile.getPeptideSeq() + " (" + str(len(self.mrmProfile.getTransitionList())) + " transitions)")        
        self.step = 100
        self.progressbar.setValue(self.step)
        self.scroll.verticalScrollBar().setSliderPosition(0)
        
    def showCurrentRow(self):
        
        self.row_index = self.peptide_list_table.currentIndex().row()
        chkVal = self.peptide_list_table.item(self.row_index,1).text()
        if chkVal == "NA":
            self.show_messageBox("정보", "해당 target의 Profiling 정보가 없습니다. 다른 Peptide를 선택하세요!")
            return
        
        self.progressbar.setValue(0)
        
        if len(self.currentWidgets) > 0:
            for w in self.currentWidgets:
                self._layout.removeWidget(w)
            self.currentWidgets.clear()            
        
        self.target_transition_list.clear()
        
        self.peptide_seq = self.peptide_list_table.item(self.row_index,0).text()
        
        self.show_msg(self.peptide_seq + " processing ...")
        
        self.mrmProfile = self.profileList[self.peptide_seq]
        
        color_code = 0
        
        min_rt = 99999999
        max_rt = 0
        
        raw_min_rt = 99999999
        raw_max_rt = 0
        
        new_transition_list = list()
        
        if self.mrmProfile.get_is_modeling() == 0:  # 만약 모델링 안된 mrmProfil이라면 모델링을 하고 했다면 그냥 뿌리가만하자 
            transitionList = self.mrmProfile.getTransitionList()
            
            transitionIndex = 0
            for transition in transitionList:
                transition_type = transition.getIonType()
                transitionIndex += 1
                # progressbar 업데이트 
                self.step = int(transitionIndex / len(transitionList) * 100)
                self.progressbar.setValue(self.step)
                self.show_msg(transition_type + " 처리 중 ...")
                
                _val = transition.checkEndoPriority()
                
                # 모델링 상관없이 new_transition_list에 등록 
                new_transition_list.append(transition)
                
                if _val == 1:
                    gen_model = transition.fit_model_reverse()
                    if raw_min_rt > transition.getRawMinRT():
                        raw_min_rt = transition.getRawMinRT()
                    if raw_max_rt < transition.getRawMaxRT():
                        raw_max_rt = transition.getRawMaxRT()
                    if gen_model == 0:
                        continue
                else:
                    gen_model = transition.fit_model()
                    if raw_min_rt > transition.getRawMinRT():
                        raw_min_rt = transition.getRawMinRT()
                    if raw_max_rt < transition.getRawMaxRT():
                        raw_max_rt = transition.getRawMaxRT()
                    if gen_model == 0:
                        continue
                    
                if len(transition.get_model_SIL()) > 0 and len(transition.get_model_Endo()) > 0:
                    model_SIL = transition.get_model_SIL()[0]
                    model_Endo = transition.get_model_Endo()[0]
                    
                    sil_model_area = model_SIL.get_model_area()
                    min_model_x, max_model_x = model_SIL.get_model_range()
                    endo_model_area = model_Endo.get_model_area_with_range(min_model_x, max_model_x)
                    ratio = round(endo_model_area / model_SIL.get_model_area(),4)
                    
                    #ratio = round(model_Endo.get_model_area() / model_SIL.get_model_area(),4)
                    rt_list_SIL = model_SIL.get_x()
                    #if max(rt_list_SIL) - min(rt_list_SIL) < 2:  # 2분 이상 인것들은 구간 측정에서 제외 
                    if min_rt > min(rt_list_SIL):
                        min_rt = min(rt_list_SIL)
                    if max_rt < max(rt_list_SIL):
                        max_rt = max(rt_list_SIL)
                    
                    #print("view:", min_rt, " ~ ", max_rt)
                    # 자동화시 아래 주석제거
                    #self.add_target_transition(transition)
    
            
            if min_rt != 99999999:
                self.mrmProfile.set_min_rt(min_rt)
                self.mrmProfile.set_max_rt(max_rt)
                self._rt_from_tf.setText(str(round(min_rt,2)))
                self._rt_to_tf.setText(str(round(max_rt,2)))
                #self.mrmProfile.set_is_modeling(1)
            else:
                self.mrmProfile.set_min_rt(raw_min_rt)
                self.mrmProfile.set_max_rt(raw_max_rt)
            
            self.mrmProfile.setModeled_TransitionList(new_transition_list)
        else:
            min_rt = self.mrmProfile.get_min_rt()
            max_rt = self.mrmProfile.get_max_rt()
            new_transition_list = self.mrmProfile.getModeled_TransitionList()
        
        transitionIndex = 1
        
        '''
        if len(new_transition_list) == 0:
            self.show_messageBox("정보", "해당 Transition의 data가 없습니다")
            return
        '''
        
        for transition in new_transition_list:
            transition.set_rt_min_for_view(self.mrmProfile.get_min_rt())
            transition.set_rt_max_for_view(self.mrmProfile.get_max_rt())
            
            # 그래프 구성 
            cWidget = CanvasWidget(self, self.mrmProfile.getPeptideSeq(), transition, transitionIndex)               
            self._layout.addWidget(cWidget)
            self.currentWidgets.append(cWidget)
            transitionIndex += 1
        
        #self.statusbar.removeWidget(self.progressbar)
        
        self.show_msg(self.mrmProfile.getPeptideSeq() + " (" + str(len(self.mrmProfile.getTransitionList())) + " transitions)")        
        self.progressbar.setValue(100)
        self.scroll.verticalScrollBar().setSliderPosition(0)
        
    def refit(self, target_transition, model_cnt_SIL, model_cnt_Endo):
        
        self.progressbar.setValue(0)
        if len(self.currentWidgets) > 0:
            for w in self.currentWidgets:
                self._layout.removeWidget(w)
            self.currentWidgets.clear()
            
        #self.threadclass.start()
        self.show_msg(target_transition.getIonType() + " processing ...")
        
        rt_from = target_transition.get_rt_min_for_view()
        rt_to   = target_transition.get_rt_max_for_view()
        
        target_transition.set_final_model_index_SIL(0)
        target_transition.set_final_model_index_Endo(0)
        
        _val = target_transition.checkEndoPriority()

        if _val == 1:
            gen_model = target_transition.fit_model_custom_reverse(model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, 1)
            if gen_model != -1:
                model_SIL = target_transition.get_model_SIL()[0]
               
        else:
            gen_model = target_transition.fit_model_custom(model_cnt_SIL, model_cnt_Endo, rt_from, rt_to, 1)
            if gen_model != -1:
                model_Endo = target_transition.get_model_Endo()[0]
            
        transitionIndex = 1
        transition_list = self.mrmProfile.getTransitionList()
        for transition in transition_list:
            
            transition.set_rt_min_for_view(self.mrmProfile.get_min_rt())
            transition.set_rt_max_for_view(self.mrmProfile.get_max_rt())
            
            # 그래프 구성 
            cWidget = CanvasWidget(self, self.mrmProfile.getPeptideSeq(), transition, transitionIndex)               
            self._layout.addWidget(cWidget)
            self.currentWidgets.append(cWidget)
            transitionIndex += 1
        
        self.progressbar.setValue(100)
        self.show_msg(target_transition.getIonType() + " fit ... Done")
        
    def change_model(self, target_transition, model_idx_SIL, model_idx_Endo):
        
        self.progressbar.setValue(0)
        if len(self.currentWidgets) > 0:
            for w in self.currentWidgets:
                self._layout.removeWidget(w)
        
        
        target_transition.set_final_model_index_SIL(model_idx_SIL)
        target_transition.set_final_model_index_Endo(model_idx_Endo)
        
        self.show_msg(target_transition.getIonType() + " processing ...")
        
        transitionIndex = 1
        transition_list = self.mrmProfile.getTransitionList()
        for transition in transition_list:
            transition.set_rt_min_for_view(self.mrmProfile.get_min_rt())
            transition.set_rt_max_for_view(self.mrmProfile.get_max_rt())
            
            # 그래프 구성 
            cWidget = CanvasWidget(self, self.mrmProfile.getPeptideSeq(), transition, transitionIndex)               
            self._layout.addWidget(cWidget)
            self.currentWidgets.append(cWidget)
            transitionIndex += 1
        
        self.progressbar.setValue(100)
        self.show_msg(target_transition.getIonType() + " change model ... Done")        
        
    def save_ratio(self):
        
        if len(self.target_transition_list) == 0:
            self.show_messageBox("Save Ratio", "선택된 transition이 없습니다!")
            return
        
        ratio_sum = 0.0
        sil_area_sum = 0.0
        endo_area_sum = 0.0
        avg_rt_sum = 0.0
        
        sil_area_str = ""
        endo_area_str = ""
        all_ratio_str = ""
        
        rep_area_SIL = 0.0
        rep_area_Endo = 0.0
        
        # transition_list에서 use_transition = 0 으로 초기화 
        transitionList = self.mrmProfile.getTransitionList()
        for transition in transitionList:
            transition.set_use_transition(0)

        plt.figure(figsize=(6, 4))
        plt.subplots_adjust(left=0.2, bottom=0.2, right=0.85, top=0.8, wspace=0.6, hspace=0.20)
        
        transitionIndex = 1
        min_x = self.mrmProfile.get_min_rt()
        max_x = self.mrmProfile.get_max_rt()
        
        img_cnt = len(self.target_transition_list) + 2
        merged_image = Image.new('RGB', (640 * img_cnt, 480), (250,250, 250))
        img_position = 0
        
        target_transition_str = ""
        
        for transition in self.target_transition_list:
            transition.set_use_transition(1)
            
            model_SIL = transition.get_model_SIL()[transition.get_final_model_index_SIL()]
            model_Endo = transition.get_model_Endo()[transition.get_final_model_index_Endo()]
            
            peak_width = model_SIL.get_peak_width()
            transition.peak_width_list.append(peak_width)
            
            
            ratio = round(model_Endo.get_model_area() / model_SIL.get_model_area(),4)
            #print(ratio)
            ratio_sum  += ratio
            sil_area_sum += model_SIL.get_model_area()
            endo_area_sum += model_Endo.get_model_area()
            
            sil_area_str += "\n" + str(round(model_SIL.get_model_area(),2))
            endo_area_str += "\n" + str(round(model_Endo.get_model_area(),2))
            all_ratio_str += "\n" + str(round(ratio, 4))
            
            # SIL이 큰 경우 
            if model_SIL.get_model_area() > model_Endo.get_model_area():
                if rep_area_SIL < model_SIL.get_model_area():
                    rep_area_SIL = model_SIL.get_model_area()
                    rep_area_Endo = model_Endo.get_model_area()
            else:
                if rep_area_Endo < model_Endo.get_model_area():
                    rep_area_SIL = model_SIL.get_model_area()
                    rep_area_Endo = model_Endo.get_model_area()                
            
            if model_SIL.get_model_center() > 0:
                avg_rt_sum += model_SIL.get_model_center()
            elif model_Endo.get_model_center() > 0:
                avg_rt_sum += model_Endo.get_model_center()
            
            target_transition_str += "\n" + transition.getIonType()
            
            
            ############### Graph 출력 #################
            profileSIL = transition.getProfileSIL()
            profileEndo = transition.getProfileEndo()
            rtList = list(profileSIL.keys())
            intensityListSIL = list(profileSIL.values())
            intensityListEndo = list(profileEndo.values())
            ax = plt.subplot(1, 1, 1)
            
            title = transition.getIonType() + "\n(" + "{:.3f}".format(float(transition.getPrecursorMz())) + "Th -> " + "{:.3f}".format(float(transition.getTransitionMz())) + "Th)\n"
            ax.set_title(title)
            ax.set_xlim(min_x, max_x)
            new_rt_list_SIL, new_intensity_list_SIL = self.check_range(rtList, intensityListSIL, min_x, max_x)
            new_rt_list_Endo, new_intensity_list_Endo = self.check_range(rtList, intensityListEndo, min_x, max_x)
            
            ax.scatter(new_rt_list_SIL, new_intensity_list_SIL, marker=".", s=10, color=self.color_list[transitionIndex])
            ax.scatter(new_rt_list_SIL, new_intensity_list_Endo, marker=".", s=10, color="blue")
            
            if len(transition.get_model_Endo()) > 0:
                model_Endo = transition.get_model_Endo()[transition.get_final_model_index_Endo()]
                model_Endo_cnt = len(transition.get_model_Endo())
                new_x_Endo = model_Endo.get_x()
                new_y_Endo = model_Endo.get_y()
                ax.plot(new_x_Endo, new_y_Endo, color="blue")                
                sse_Endo = model_Endo.get_sse()
                model_area_Endo = model_Endo.get_model_area()
                
                
        
            #print("Model SIL:", self.transition.get_model_SIL(), ":", self.transition.get_model_Endo())
            if len(transition.get_model_SIL()) > 0:
                model_SIL = transition.get_model_SIL()[transition.get_final_model_index_SIL()]
                model_SIL_cnt = len(transition.get_model_SIL())
                new_x_SIL = model_SIL.get_x()
                new_y_SIL = model_SIL.get_y()
                ax.plot(new_x_SIL, new_y_SIL, color=self.color_list[transitionIndex])                
                sse_SIL = model_SIL.get_sse()
                model_area_SIL = model_SIL.get_model_area()
                #ratio = round(self.model_area_Endo / self.model_area_SIL,4)               
                #x_label = "ratio:" + str(ratio) + "\nSSE(SIL):" + str(sse_SIL) + "\nSSE(Endo):" + str(sse_Endo)
                #x_label = "SSE(SIL):" + str(sse_SIL) + "\nSSE(Endo):" + str(sse_Endo)
            
            ####### 서브 plot #######
            if model_SIL.get_model_area() > model_Endo.get_model_area():
                axins = inset_axes(ax, width='30%', height='30%', loc='upper right', borderpad=2.0)
                axins.scatter(new_rt_list_Endo, new_intensity_list_Endo, marker=".", s=10, color="blue")
                if max(new_intensity_list_Endo) > 0:
                    axins.plot(new_x_Endo, new_y_Endo, label='Graph 2', color='blue')
                    
                #axins.set_xlabel('')
                #axins.set_ylabel('', color='r')
                axins.tick_params(axis='y', labelcolor='r')
            else:
                axins = inset_axes(ax, width='30%', height='30%', loc='upper right', borderpad=2.0)
                axins.scatter(new_rt_list_SIL, new_intensity_list_SIL, marker=".", s=10, color=self.color_list[transitionIndex])
                if max(new_intensity_list_SIL) > 0:
                    axins.plot(new_x_SIL, new_y_SIL, label='Graph 2', color=self.color_list[transitionIndex])
                #axins.set_xlabel('')
                #axins.set_ylabel('', color='r')
                axins.tick_params(axis='y', labelcolor='r')
            
            
            x_label = "RT"
            ax.set_xlabel(x_label)
            full_path = self.write_img(self.mrmProfile.getPeptideSeq(), transition.getIonType())            
            transition.set_img_path(full_path)
            plt.clf()
            transitionIndex += 1
            
            ######### 이미지 병합 ###########
            image = Image.open(full_path)
            merged_image.paste(image, (img_position, 0))
            img_position += 640
         
        ####################### 모든 transition을 하나로 합친 거 #####################
        transitionIndex = 1
        plt.figure(figsize=(12, 4))
        plt.subplots_adjust(top=0.8)
        entire_sil_ax = plt.subplot(1, 2, 1)
        entire_endo_ax = plt.subplot(1, 2, 2)
        
        entire_sil_ax.set_title(self.mrmProfile.getPeptideSeq() + "\n(SIL)\n")
        entire_endo_ax.set_title(self.mrmProfile.getPeptideSeq() + "\n(Endo)\n")
        
        entire_sil_ax.set_xlabel("RT")
        entire_endo_ax.set_xlabel("RT")
        
        entire_sil_ax.set_xlim(min_x, max_x)
        entire_endo_ax.set_xlim(min_x, max_x)
        
        for transition in self.target_transition_list:
            transition.set_use_transition(1)
            
            ############### Graph 출력 #################
            profileSIL = transition.getProfileSIL()
            profileEndo = transition.getProfileEndo()
            rtList = list(profileSIL.keys())
            intensityListSIL = list(profileSIL.values())
            intensityListEndo = list(profileEndo.values())
            
            new_rt_list_SIL, new_intensity_list_SIL = self.check_range(rtList, intensityListSIL, min_x, max_x)
            new_rt_list_Endo, new_intensity_list_Endo = self.check_range(rtList, intensityListEndo, min_x, max_x)
            
            entire_sil_ax.scatter(new_rt_list_SIL, new_intensity_list_SIL, marker=".", s=3, color=self.color_list[transitionIndex])
            entire_endo_ax.scatter(new_rt_list_SIL, new_intensity_list_Endo, marker=".", s=3, color=self.color_list[transitionIndex])
            
            if len(transition.get_model_Endo()) > 0:
                model_Endo = transition.get_model_Endo()[transition.get_final_model_index_Endo()]
                model_Endo_cnt = len(transition.get_model_Endo())
                new_x_Endo = model_Endo.get_x()
                new_y_Endo = model_Endo.get_y()                
                entire_endo_ax.plot(new_x_Endo, new_y_Endo, color=self.color_list[transitionIndex])                
        
            #print("Model SIL:", self.transition.get_model_SIL(), ":", self.transition.get_model_Endo())
            if len(transition.get_model_SIL()) > 0:
                model_SIL = transition.get_model_SIL()[transition.get_final_model_index_SIL()]
                model_SIL_cnt = len(transition.get_model_SIL())
                new_x_SIL = model_SIL.get_x()
                new_y_SIL = model_SIL.get_y()
                entire_sil_ax.plot(new_x_SIL, new_y_SIL, color=self.color_list[transitionIndex])                
                
            transitionIndex += 1
            
        ######### 이미지 병합 ###########
        full_path = self.entire_write_img(self.mrmProfile.getPeptideSeq())    
        image = Image.open(full_path)
        merged_image.paste(image, (img_position, 0))
        img_position += 640            
        
        ########### 각 transition별로 생성된 이미지를 병합 #########
        merged_image_path = self.img_dir_name + "/" + self.mrmProfile.getPeptideSeq() + ".png"
        
        
        try:
            merged_image.save(merged_image_path)
            self.mrmProfile.setImgPath(merged_image_path)
        except OSError as e:    
                print(e)
        
        
        #avg_ratio = round(ratio_sum / len(self.target_transition_list),4)
        avg_sil_area  = round(sil_area_sum / len(self.target_transition_list),4)
        avg_endo_area  = round(endo_area_sum / len(self.target_transition_list),4)
        avg_ratio = avg_endo_area / avg_sil_area
        avg_ratio = round(avg_ratio, 4)
        avg_rt = round(avg_rt_sum / len(self.target_transition_list),4)
        
        self.mrmProfile.set_rep_area_SIL(rep_area_SIL)
        self.mrmProfile.set_rep_area_Endo(rep_area_Endo)
        
        self.mrmProfile.set_avg_ratio(avg_ratio)
        self.mrmProfile.set_avg_SIL_Area(avg_sil_area)
        self.mrmProfile.set_avg_Endo_Area(avg_endo_area)
        if len(sil_area_str) > 0:
            self.mrmProfile.set_avg_SIL_Area_Str(sil_area_str[1:])
        else:
            self.mrmProfile.set_avg_SIL_Area_Str("0.0")
            
        if len(endo_area_str) > 0:
            self.mrmProfile.set_avg_Endo_Area_Str(endo_area_str[1:])
        else:
            self.mrmProfile.set_avg_Endo_Area_Str("0.0")
            
        if len(all_ratio_str) > 0:
            self.mrmProfile.set_all_ratio(all_ratio_str[1:])
        else:
            self.mrmProfile.set_all_ratio("0.0")
            
        self.mrmProfile.set_avg_rt(avg_rt)
        
        if len(target_transition_str) > 0:
            self.mrmProfile.set_target_transition_str(target_transition_str[1:])
        else:
            self.mrmProfile.set_target_transition_str("")
        #print("avg ratio:", avg_ratio)
        self.peptide_list_table.setItem(self.row_index, 2, QTableWidgetItem(str(avg_ratio)))
        self.peptide_list_table.item(self.row_index, 0).setBackground(QtGui.QColor(153,204,255))
        self.peptide_list_table.item(self.row_index, 1).setBackground(QtGui.QColor(153,204,255))
        self.peptide_list_table.item(self.row_index, 2).setBackground(QtGui.QColor(153,204,255))
        self.show_msg("저장되었습니다")
        
        self.mrmProfile.set_is_modeling(1)
        ########### Threading으로 Excel에 저장 ###########
        
        self.saveCnt += 1
        
        if self.saveCnt >= 5:
            try:
                output_file_name = self.input_file_name[0: self.input_file_name.rfind(".")] + "_Excel.xlsx"
                print("saving Excel ... ", output_file_name)
                self.make_excel(output_file_name)
            except xlsxwriter.exceptions.FileCreateError as e:
                self.show_messageBox("경고", "결과를 저장할 수 없습니다. Excel 파일이 열려 있는지 확인하세요")
                #self.show_msg("결과를 저장할 수 없습니다. Excel 파일일 오픈되어 있는지 확인하세요")
            self.saveCnt = 0
            
        
    def write_img(self, peptide_seq, transitionType):        
        
        self.img_dir_name = self.input_file_name[0: self.input_file_name.rfind(".")]
        img_dir = self.img_dir_name + "/"
        try:
            if not os.path.exists(img_dir):
                os.makedirs(img_dir)
        except OSError:
                print("Failed to create the image directory")
        img_name = peptide_seq + "_" + transitionType + ".png"
        full_path = img_dir + "/" + img_name
                
        #print(full_path, " saving ...")
        try:
            plt.savefig(full_path)
        except OSError as e:    
                print(e)
        
        return full_path
    
    def entire_write_img(self, peptide_seq):        
        
        self.img_dir_name = self.input_file_name[0: self.input_file_name.rfind(".")]
        img_dir = self.img_dir_name + "/"
        try:
            if not os.path.exists(img_dir):
                os.makedirs(img_dir)
        except OSError:
                print("Failed to create the image directory")
        img_name = peptide_seq + "_All.png"
        full_path = img_dir + "/" + img_name
                
        #print(full_path, " saving ...")
        try:
            plt.savefig(full_path)
        except OSError as e:    
                print(e)
        
        return full_path    
    
    def check_range(self, rtList, intensityList, min_x, max_x):
        new_rt_list = list()
        new_intensity_list = list()
        
        for i in range(0, len(rtList)):
            rt = rtList[i]
            if rt < min_x or rt > max_x:
                continue
            new_rt_list.append(rt)
            new_intensity_list.append(intensityList[i])
        
        return new_rt_list, new_intensity_list
        
    def show_messageBox(self, title, msg):
        QMessageBox.about(self, title, msg)
        
    def make_excel(self, output_file_name):
        color_list = [ "red", "green", "orange", "purple", "brown", "pink", "gray", "#808000", "#00FFFF", "#8B0000", "#006400", "#FF8C00", "#EE82EE", "#8B4513", "#FF00FF", "#008B8B", "#8B008B" ]
        
        
        wb = xlsxwriter.Workbook(output_file_name)
        worksheet = wb.add_worksheet()
        cell_format = wb.add_format()
        cell_format.set_center_across()
        cell_format.set_align('vcenter')
        cell_format.set_text_wrap()
        
        font_format_01 = wb.add_format();
        font_format_01.set_center_across()
        font_format_01.set_align("vcenter")
        font_format_01.set_text_wrap()
        font_format_01.set_font_color(color_list[0])
        
        font_format_02 = wb.add_format();
        font_format_02.set_center_across()
        font_format_02.set_align("vcenter")
        font_format_02.set_text_wrap()
        font_format_02.set_font_color(color_list[1])
        
        font_format_03 = wb.add_format();
        font_format_03.set_center_across()
        font_format_03.set_align("vcenter")
        font_format_03.set_text_wrap()
        font_format_03.set_font_color(color_list[2])
        
        font_format_04 = wb.add_format();
        font_format_04.set_center_across()
        font_format_04.set_align("vcenter")
        font_format_04.set_text_wrap()
        font_format_04.set_font_color(color_list[3])
        
        font_format_05 = wb.add_format();
        font_format_05.set_center_across()
        font_format_05.set_align("vcenter")
        font_format_05.set_text_wrap()
        font_format_05.set_font_color(color_list[4])
        
        font_format_06 = wb.add_format();
        font_format_06.set_center_across()
        font_format_06.set_align("vcenter")
        font_format_06.set_text_wrap()
        font_format_06.set_font_color(color_list[5])
        
        font_format_07 = wb.add_format();
        font_format_07.set_center_across()
        font_format_07.set_align("vcenter")
        font_format_07.set_text_wrap()
        font_format_07.set_font_color(color_list[6])
        
        font_format_08 = wb.add_format();
        font_format_08.set_center_across()
        font_format_08.set_align("vcenter")
        font_format_08.set_text_wrap()
        font_format_08.set_font_color(color_list[7])
        
        font_format_09 = wb.add_format();
        font_format_09.set_center_across()
        font_format_09.set_align("vcenter")
        font_format_09.set_text_wrap()
        font_format_09.set_font_color(color_list[8])
        
        font_format_10 = wb.add_format();
        font_format_10.set_center_across()
        font_format_10.set_align("vcenter")
        font_format_10.set_text_wrap()
        font_format_10.set_font_color(color_list[9])
        
        font_format_11 = wb.add_format();
        font_format_11.set_center_across()
        font_format_11.set_align("vcenter")
        font_format_11.set_text_wrap()
        font_format_11.set_font_color(color_list[10])
        
        font_format_12 = wb.add_format();
        font_format_12.set_center_across()
        font_format_12.set_align("vcenter")
        font_format_12.set_text_wrap()
        font_format_12.set_font_color(color_list[11])
        
        font_format_13 = wb.add_format();
        font_format_13.set_center_across()
        font_format_13.set_align("vcenter")
        font_format_13.set_text_wrap()
        font_format_13.set_font_color(color_list[12])
        
        font_format_14 = wb.add_format();
        font_format_14.set_center_across()
        font_format_14.set_align("vcenter")
        font_format_14.set_text_wrap()
        font_format_14.set_font_color(color_list[13])
        
        font_format_15 = wb.add_format();
        font_format_15.set_center_across()
        font_format_15.set_align("vcenter")
        font_format_15.set_text_wrap()
        font_format_15.set_font_color(color_list[14])
        
        font_format_list = list()
        font_format_list.append(font_format_01)
        font_format_list.append(font_format_02)
        font_format_list.append(font_format_03)
        font_format_list.append(font_format_04)
        font_format_list.append(font_format_05)
        font_format_list.append(font_format_06)
        font_format_list.append(font_format_07)
        font_format_list.append(font_format_08)
        font_format_list.append(font_format_09)
        font_format_list.append(font_format_10)
        font_format_list.append(font_format_11)
        font_format_list.append(font_format_12)
        font_format_list.append(font_format_13)
        font_format_list.append(font_format_14)
        font_format_list.append(font_format_15)
        
        
        worksheet.set_column("A:A", 50)
        worksheet.set_column("B:B", 20)
        worksheet.set_column("C:C", 20)
        worksheet.set_column("D:D", 20)
        worksheet.set_column("E:E", 20)
        worksheet.set_column("F:F", 20)
        worksheet.set_column("G:G", 20)
        worksheet.set_column("H:H", 20)
        worksheet.set_column("I:I", 20)
        worksheet.set_column("J:J", 350)
        
        
        worksheet.write("A1", "Peptide", cell_format)
        worksheet.write("B1", "RT", cell_format)  
        worksheet.write("C1", "Precursor m/z(SIL)", cell_format)  
        worksheet.write("D1", "Transition", cell_format)
        worksheet.write("E1", "Area(SIL)", cell_format)
        worksheet.write("F1", "Area(Endo)", cell_format)
        worksheet.write("G1", "Ratio", cell_format)
        worksheet.write("H1", "Rep Area(SIL)", cell_format)
        worksheet.write("I1", "Rep Area(Endo)", cell_format)
        worksheet.write("J1", "Profile", cell_format)
        
        row_index = 2
        cell_names = ("A","B","C","D","E","F","G","H","I","J")
        #cell_names = ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF","AG", "AH", "AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV")
        
        peptideList = self.profileList.keys()
        
        i = 0
        
        sim_data = list()
        sim_data.append("Peptide\tm/z Start\tm/z End\tRT\tpeak_width\n")
        
        for peptide_seq in peptideList:
            
            mrmProfile = self.profileList[peptide_seq]
            transitionList = mrmProfile.getTransitionList()
            
            if len(transitionList) == 0:
                continue
            
            if mrmProfile.img_path == None:
                continue
            
            # SIM 수집 
            precursor_mz = mrmProfile.getPrecursorMz()
            from_mz = float(precursor_mz) - 4.01
            from_mz = round(from_mz, 4)
            to_mz   = float(precursor_mz) + 1.0
            to_mz = round(to_mz, 4)
            peak_rt = round(mrmProfile.get_avg_rt(), 2)
            peak_width = mrmProfile.get_peak_width()
            sim_info = mrmProfile.getPeptideSeq() + "\t" + str(from_mz) + "\t" + str(to_mz) + "\t" + str(peak_rt) + "\t" + str(peak_width) + "\n"
            sim_data.append(sim_info)
            
            worksheet.set_row ((i+1), 153)
            
            cell_index = cell_names[0] + str(row_index)
            worksheet.write(cell_index, mrmProfile.getPeptideSeq(), cell_format)
            
            cell_index = cell_names[1] + str(row_index)
            worksheet.write(cell_index, mrmProfile.get_avg_rt(), cell_format)
            
            
            cell_index = cell_names[2] + str(row_index)
            worksheet.write(cell_index, mrmProfile.getPrecursorMz(), cell_format)
            
            cell_index = cell_names[3] + str(row_index)
            worksheet.write(cell_index, mrmProfile.get_target_transition_str(), cell_format)
            
            cell_index = cell_names[4] + str(row_index)
            worksheet.write(cell_index, mrmProfile.get_avg_SIL_Area_Str(), cell_format)
            
            cell_index = cell_names[5] + str(row_index)
            worksheet.write(cell_index, mrmProfile.get_avg_Endo_Area_Str(), cell_format)
            
            #ratio = mrmProfile.get_avg_Endo_Area() / mrmProfile.get_avg_SIL_Area()
            #ratio = round(ratio, 4)
            
            cell_index = cell_names[6] + str(row_index)
            #worksheet.write(cell_index, str(ratio), cell_format)
            worksheet.write(cell_index, mrmProfile.get_all_ratio(), cell_format)
            
            cell_index = cell_names[7] + str(row_index)
            worksheet.write(cell_index, mrmProfile.get_rep_area_SIL(), cell_format)
            
            cell_index = cell_names[8] + str(row_index)
            worksheet.write(cell_index, mrmProfile.get_rep_area_Endo(), cell_format)
            
            cell_index = cell_names[9] + str(row_index)
            worksheet.insert_image(cell_index, mrmProfile.getImgPath(), {'x_scale':0.4, 'y_scale':0.4, 'x_offset':10, 'y_offset':10, 'object_position': 1})

            row_index += 1
            
            i += 1
        
        sim_data_file = output_file_name[0:output_file_name.rfind(".")] + "_SIMList.txt"
        
        with open(sim_data_file, 'w') as file:
            file.writelines(sim_data)
        
        wb.close()
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    sys.exit(app.exec_())
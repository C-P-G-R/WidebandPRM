# WidebandPRM

Program to find the area after fitting for data analyzed by WidebandPRM

Installation

1. download anaconda
2. conda create -n WidebandPRM python=3.7.9
3. conda activate WidebandPRM
4. install packages
5. 
     conda install numpy
   
     conda install pandas
   
     conda install spyder
   
     conda install scipy
   
     conda install matplotlib
   
     pip install lmfit
   
     pip install peakutils
   
     pip install xlsxwriter
   
     pip install image
   
   
7. run
 
   python WidebandPRM.py

Notes

The input file uses mzML.
- Use the --32 option to generate mzML
- You should have the following files in the same location as your mzML- 
- target1.csv, target1_trigger.csv (for +2)
- target2.csv, target2_trigger.csv (for +3)


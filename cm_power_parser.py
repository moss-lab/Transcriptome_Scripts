"""
cm_power_parser.py
By Collin A. O'Leary, Moss Lab, ISU
Usage: python cm_power_parser.py >> output_file_name.txt
run in directory with power files and it will extract the covarying basepairs and bin them based on the power threshold
"""

import sys
import os
import glob


CWD = os.getcwd()         
for x in os.walk(CWD):  #Walks through working directory
    pwr = glob.glob("*.power")  #Makes a list with all power files
print(f"Name\tBPs\tavg_substitutions\texpected BPs\tSTD_DEV\tObserved BPs\t# of BP power 0-0.1\t# of BP power 0.1-0.25\t# of BP power >=0.25\t0-0.1 BP info\t0.1-0.25 BP info\t>=0.25 BP info") #prints the file header
pwr_data = []    
for obj in pwr:             #parse throughs the list of power files
    with open(obj) as f:    #opens individual power file
        lines = f.readlines()
        name = obj.replace("_1.power","")
        low_cv = []
        mid_cv = []
        high_cv = []
        for x in range(0,len(lines)):   #loop that parses through lines of power file
            line = lines[x].split()
            if len(line) <= 1:      #if loop which grabs different bits of data, depending on what line is being parsed
                pass
            elif line[1] == "BPAIRS" and line[2] == "expected":
                ex_bp = float(line[5])
                std_dev = float(line[7])
            elif line[1] == "BPAIRS" and line[2] == "observed":
                obv_bp = float(line[5])
            elif line[1] == "BPAIRS":
                bps = float(line[2])
            elif line[1] == "avg":
                subs = float(line[5])
            elif line[0] == "*":            #this sub loop seperates out the data based on the degree of power
                if float(line[4]) >= 0.25:
                    high_cv.append((int(line[1]),int(line[2]),float(line[4])))
                elif float(line[4]) >= 0.1:
                    mid_cv.append((int(line[1]),int(line[2]),float(line[4])))
                else:
                    low_cv.append((int(line[1]),int(line[2]),float(line[4])))
            else:
                pass
    low_num = len(low_cv)
    mid_num = len(mid_cv)
    high_num = len(high_cv)
    print(f"{name}\t{bps}\t{subs}\t{ex_bp}\t{std_dev}\t{obv_bp}\t{low_num}\t{mid_num}\t{high_num}\t{low_cv}\t{mid_cv}\t{high_cv}")   #prints resulting data
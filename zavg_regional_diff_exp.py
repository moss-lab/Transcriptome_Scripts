#zavg_regional_diff_exp
"""
By Collin A. O'Leary
Usage: python zavg_regional_diff_exp.py <output from regional_zavg.py> <list of gene names to parse data to>

"""

import sys
import os
import glob
import statistics
from statistics import mean
from datetime import datetime 

filename1 = sys.argv[1]	        #This is transcriptome zavg metrics output
filename2 = sys.argv[2]         #This is a list of sub genes to pull data out for, single column name only

tran_met_dic = {}                   #A dictionary for transcriptome metrics
with open(filename1 , 'r') as tm:   #Opens the zavg regional results file
    tran_data = tm.readlines()
    for x in range(1,len(tran_data)):   #A loop that parses the zavg regional data and splits it apart
        z = tran_data[x].split()
        y = z[0].split("_")
        for obj in y:                   #A loop to define the ENST corresponding to the data
            if "ENST" in obj and "." in obj:
                name = obj
                
            else:
                pass
        tran_met_dic[name] = z      #Puts the enst id and data into a dictionary
    
with open(filename2 , 'r') as st:   #Opens the file with gene names to parse data to
    sub_trans = st.readlines()
    

out_tran = []

for line in sub_trans:          #Iterates through the list of gene names
    t = line.split()
    if t[0] in tran_met_dic:        #If the gene name is in the dictionary, it pulls out the data to a sub list
        tran = tran_met_dic[t[0]]
        out_tran.append(tran)
    else:
        pass



current_date = datetime.now().date()

out_name = filename2.replace('.txt','')
outfile1 = f"Regional_Zavg_of_{out_name}_{current_date}.txt"    #Grabs the data and time to time stamp the outfile

with open(outfile1, "w") as tm: #Opens the out file
    for obj in out_tran:
        
        tm.write(f"{obj[0]}\t {obj[1]}\t{obj[2]}\t{obj[3]}\t{obj[4]}\t{obj[5]}\t{obj[6]}\n")    #Writes the data from the sub list of data to an output file
        
            

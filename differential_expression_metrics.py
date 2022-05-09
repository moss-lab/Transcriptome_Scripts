"""
differential_expression_metrics
By Collin A. O'Leary
Usage: python differential_expression_metrics.py <transcriptome_metrics output> <list of genes to parse> 
run this with the full transcriptome metrics output from ScanFold. It also needs a sub list of genes to check for and cross reference
"""

import sys
import os
import glob
import statistics
from statistics import mean
from datetime import datetime 

filename1 = sys.argv[1]	        #This is transcriptome metrics output
filename2 = sys.argv[2]         #This is a list of sub genes to pull data out for

tran_met_dic = {}           #creates a dictionary for transcriptome metrics
with open(filename1 , 'r') as tm:   #opens transcriptome metric output, makes a dictionary with the enst ID as a key and attaches the data as a value
    tran_data = tm.readlines()
    for x in range(1,len(tran_data)):
        z = tran_data[x].split()
        y = z[0].split("_")
        for obj in y:
            if "ENST" in obj and "." in obj:
                name = obj
            else:
                pass
        tran_met_dic[name] = z
    
with open(filename2 , 'r') as st:   #opens the transcript ID file and makes it a list
    sub_trans = st.readlines()
    

out_tran = []

for line in sub_trans:  #Parses the transcript list to pull out data from the dictionary, if present
    t = line.split()
    if t[0] in tran_met_dic:
        tran = tran_met_dic[t[0]]
        out_tran.append(tran)
    else:
        pass
        



current_date = datetime.now().date()

out_name = filename2.replace('.txt','')
outfile1 = f"Transcriptome_metrics_of_{out_name}_{current_date}.txt"    #grabs date and time for accurate out file tracking

with open(outfile1, "w") as tm:     #writes the parsed output data to a file
    for obj in out_tran:
        tm.write(f"{obj[0]}\t {obj[1]}\t{obj[2]}\t{obj[3]}\t{obj[4]}\t{obj[5]}\t{obj[6]}\t{obj[7]}\t{obj[8]}\t{obj[9]}\t{obj[10]}\n")

"""
regional_zavg
By Collin A. O'Leary, Moss Lab, ISU
Usage: python regional_zavg.py <GFF3 file>
Compares a GFF3 file to all available ScanFold zavg wig files from ScanFold output directories and parses them to appropriate regions
"""

import sys
import os
import glob
import statistics
from statistics import mean
from datetime import datetime 

filename1 = sys.argv[1]	        #GFF3 with CDS and UTR info
region_dict = {}
with open(filename1 , 'r') as gff3: #Opens the GFF3 file to parse the data
    lines = gff3.readlines()
    for x in range(0, len(lines)):      #this loop parses each line of the file looking for unique identifies which allow extraction of regional info to a dictionary
        if lines[x].startswith("##"):
            y = lines[x].split()
            enst = y[1]
            z = lines[x+1].split()
            v = z[8].split(";")
            
            if "UTR5" in z[8]:
                k = z[8].split(";")
                for obj in k:
                    if obj.startswith("UTR5"):
                        g = obj.split(":")
                        UTR5 = g[1]
                    else:
                        pass
            else:
                UTR5 = "0-0"

            if "CDS" in z[8]:
                k = z[8].split(";")
                for obj in k:
                    if obj.startswith("CDS:"):
                        g = obj.split(":")
                        
                        CDS = g[1]
                    else:
                        pass
            else:
                CDS = "0-0"
        
            if "UTR3" in z[8]:
                k = z[8].split(";")
                for obj in k:
                    if obj.startswith("UTR3"):
                        g = obj.split(":")
                        UTR3 = g[1]
                    else:
                        pass
            else:
                UTR3 = "0-0"
            region_dict[enst] = (UTR5,CDS,UTR3)
        else:
            pass


CWD = os.getcwd()
dirs = []           #List of all directories with (ENST, root) format, skip position 0

for x in os.walk(CWD):      #Walks through directory to build a list of paths
    
    y = x[0]
    z = y.replace(CWD+"/","")    #Need to fix this replace command depending on whether its run on windows or linux
    u = z.split("_")
    for obj in u:
        if "ENST" in obj and "." in obj:
            dirs.append((obj,x[0],z))
        else:
            pass

enst_data = []
corrupted_dirs = []

def avg_zs(list):       #Used for all list averaging
    if len(list) == 0:
        return ''
    else:
        return statistics.mean(list)

for x in range(0,len(dirs)):        #parses through directory list, grabs the Zavg wig file from ScanFold and extracts z-score data
    
    os.chdir(dirs[x][1])
    zavgfile = glob.glob("*_Zavg_metrics*")
    if dirs[x][0] in region_dict:
        regions = region_dict[dirs[x][0]]
    else:
        corrupted_dirs.append(dirs[x])
        continue
    
    if len(zavgfile) == 0:
        corrupted_dirs.append(dirs[x])
        continue
    twd = os.getcwd()
    #print(twd)
    with open(zavgfile[0], "r") as of:
        lines = of.readlines()
        zs_5utr = []
        zs_3utr = []
        zs_cds = []
        utr5 = regions[0].split("-")
        cds = regions[1].split("-")
        utr3 = regions[2].split("-")
        #print(utr5,cds,utr3)
        
        if utr5[0] == "0" and utr5[1] == "0":
            pass
        else:
            for y in range(int(utr5[0]),int(utr5[1])+1):
                zs_5utr.append(float(lines[y].replace('\n','')))
        #print(zs_5utr)
        
        if cds[0] == "0" and cds[1] == "0":
            pass
        else:
            for y in range(int(cds[0]),int(cds[1])+1):
                zs_cds.append(float(lines[y].replace('\n','')))
        #print(zs_cds)
        
        if utr3[0] == "0" and utr3[1] == "0":
            pass
        else:
            for y in range(int(utr3[0]),int(utr3[1])+1):
                zs_3utr.append(float(lines[y].replace('\n','')))
        #print(zs_3utr)
        #Averages the data for each list
        avg5 = avg_zs(zs_5utr)
        avgc = avg_zs(zs_cds)
        avg3 = avg_zs(zs_3utr)
        enst_data.append((dirs[x][2],zs_5utr,zs_cds,zs_3utr,avg5,avgc,avg3,regions[0],regions[1],regions[2]))
        
#print(enst_data)

os.chdir(CWD)
current_date = datetime.now().date()
current_time = datetime.now().time()

outfile1 = f"Transcriptome_regional_avg_zs_{current_date}.txt"      #date and time are grabbed for accurate output tracking

with open(outfile1, "w") as tm:     #The data is written to an output file
    tm.write("ENST\t5UTR_region\t5UTR_Avg_ZS\tCDS_region\tCDS_Avg_ZS\t3UTR_region\t3UTR_Avg_ZS\t5UTR_values\tCDS_values\t3UTR_values\n")
    for obj in enst_data:
        tm.write(f"{obj[0]}\t {obj[7]}\t{obj[4]}\t{obj[8]}\t{obj[5]}\t{obj[9]}\t{obj[6]}\n")
outfile4 = f"Transcriptome_Metrics_corrupted_directories_{current_date}.txt"
with open(outfile4, "w") as cd:         #directories that could not be computed are written here, for quality control
    for obj in corrupted_dirs:
        cd.write(obj[0]+" " +obj[1]+"\n")

"""
regional_metrics
By Collin A. O'Leary, Moss Lab, ISU
Usage: python regional_metrics.py <GFF3 file>
 
This script takes a GFF3 file and extracts 5'UTR, CDS, and 3'UTR coordinates for each present gene
It puts location data into a dictionary
Metrics are then pulled from the ScanFold final partners file and cross referenced to regions to gain insights on regional metrics
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

for x in os.walk(CWD): #Walks through directory to build a list of paths
    
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

def avg_data(list):         #Used for all list averaging
    if len(list) == 0:
        return ''
    else:
        return statistics.mean(list)


for x in range(0,len(dirs)):    #parses through directory list, grabs the final partner file from ScanFold and extracts MFE, ED, and z-score data
    
    os.chdir(dirs[x][1])
    fpf = glob.glob("*FinalPartners*")
    if dirs[x][0] in region_dict:
        regions = region_dict[dirs[x][0]]
    else:
        corrupted_dirs.append(dirs[x])
        continue
    
    if len(fpf) == 0:
        corrupted_dirs.append(dirs[x])
        continue
    twd = os.getcwd()
    #print(twd)
    with open(fpf[0], "r") as fp:
        line = fp.readlines()[1::]
        zs_5utr = []
        zs_3utr = []
        zs_cds = []
        mfe_5utr = []
        mfe_3utr = []
        mfe_cds = []
        ed_5utr = []
        ed_3utr = []
        ed_cds = []
        utr5 = regions[0].split("-")
        cds = regions[1].split("-")
        utr3 = regions[2].split("-")
        
        
        if utr5[0] == "0" and utr5[1] == "0":
            pass
        else:
            for y in range(int(utr5[0]),int(utr5[1])+1):
                if int(y) >= int(len(line)):
                    pass
                else:
                    lines = line[y].split()
                    zs_5utr.append(float(lines[4].replace('\n','')))
        #print(zs_5utr)
            for y in range(int(utr5[0]),int(utr5[1])+1):
                if int(y) >= int(len(line)):
                    pass
                else:
                    lines = line[y].split()
                    mfe_5utr.append(float(lines[3].replace('\n','')))
            for y in range(int(utr5[0]),int(utr5[1])+1):
                if int(y) >= int(len(line)):
                    pass
                else:
                    lines = line[y].split()
                    ed_5utr.append(float(lines[5].replace('\n','')))
        
        if cds[0] == "0" and cds[1] == "0":
            pass
        else:
            for y in range(int(cds[0]),int(cds[1])+1):
                if int(y) >= int(len(line)):
                    pass
                    
                else:
                    lines = line[y].split()
                    zs_cds.append(float(lines[4].replace('\n','')))
        #print(zs_cds)
            for y in range(int(cds[0]),int(cds[1])+1):
                if int(y) >= int(len(line)):
                    pass
                    
                else:
                    lines = line[y].split()
                    mfe_cds.append(float(lines[3].replace('\n','')))
            for y in range(int(cds[0]),int(cds[1])+1):
                if int(y) >= int(len(line)):
                    pass
                    
                else:
                    lines = line[y].split()
                    ed_cds.append(float(lines[5].replace('\n','')))
                
        if utr3[0] == "0" and utr3[1] == "0":
            pass
        else:
            #print(dirs[x][1])
            #print(utr3[0],utr3[1])
            #print(len(line))
            for y in range(int(utr3[0]),int(utr3[1])+1):
                if int(y) >= int(len(line)):
                    pass
                    
                else:
                    lines = line[y].split()
                    zs_3utr.append(float(lines[4].replace('\n','')))
        #print(zs_3utr)
            for y in range(int(utr3[0]),int(utr3[1])+1):
                if int(y) >= int(len(line)):
                    pass
                    
                else:
                    lines = line[y].split()
                    mfe_3utr.append(float(lines[3].replace('\n','')))
            for y in range(int(utr3[0]),int(utr3[1])+1):
                if int(y) >= int(len(line)):
                    pass
                    
                else:
                    lines = line[y].split()
                    ed_3utr.append(float(lines[5].replace('\n','')))
        #after appropriate data is put into lists, the average of the lists are calculated below
        
        avg5 = avg_data(zs_5utr)
        avgc = avg_data(zs_cds)
        avg3 = avg_data(zs_3utr)
        avg5m = avg_data(mfe_5utr)
        avgcm = avg_data(mfe_cds)
        avg3m = avg_data(mfe_3utr)
        avg5e = avg_data(ed_5utr)
        avgce = avg_data(ed_cds)
        avg3e = avg_data(ed_3utr)
        
        zs = zs_3utr+zs_5utr+zs_cds
        ed = ed_3utr+ed_5utr+ed_cds
        mfe = mfe_3utr+mfe_5utr+mfe_cds
        
        zs_tot = avg_data(zs)
        ed_tot = avg_data(ed)
        mfe_tot = avg_data(mfe)
        
        
        
        enst_data.append((dirs[x][2],zs_tot,ed_tot,mfe_tot,avg5,avgc,avg3,regions[0],regions[1],regions[2],mfe_5utr,mfe_cds,mfe_3utr,avg5m,avgcm,avg3m,ed_5utr,ed_cds,ed_3utr,avg5e,avgce,avg3e)) #list of all compiled data
        
#print(enst_data)

os.chdir(CWD)
current_date = datetime.now().date()
current_time = datetime.now().time()

outfile1 = f"Transcriptome_regional_avg_zs_mfe_ed_{current_date}.txt"   #date and time are grabbed for accurate output tracking
#ADD MFE AND ED HEADERS AND DATA PLACES
with open(outfile1, "w") as tm:     #The data is written to an output file
    tm.write("ENST\t5UTR_region\t5UTR_Avg_ZS\t5UTR_Avg_MFE\t5UTR_Avg_ED\tCDS_region\tCDS_Avg_ZS\tCDS_Avg_MFE\tCDS_Avg_ED\t3UTR_region\t3UTR_Avg_ZS\t3UTR_Avg_MFE\t3UTR_Avg_ED\tAvg. Tot. ZS\t/avg. Tot. ED\tAvg. Tot. MFE\n")
    for obj in enst_data:
        tm.write(f"{obj[0]}\t {obj[7]}\t{obj[4]}\t{obj[13]}\t{obj[19]}\t{obj[8]}\t{obj[5]}\t{obj[14]}\t{obj[20]}\t{obj[9]}\t{obj[6]}\t{obj[15]}\t{obj[21]}\t{obj[1]}\t{obj[2]}\t{obj[3]}\n")
outfile4 = f"Transcriptome_Metrics_corrupted_directories_{current_date}.txt"
with open(outfile4, "w") as cd:     #directories that could not be computed are written here, for quality control
    for obj in corrupted_dirs:
        cd.write(obj[0]+" " +obj[1]+"\n")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

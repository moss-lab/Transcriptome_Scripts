"""
Transcriptome_Metrics
By Collin A. O'Leary, Moss Lab, ISU
Usage: python transcriptome_metrics.py
run this in a directory which contains subdirectories of ScanFold output to get metrics on the datasets also checks whether output is present with complete data
"""

import sys
import os
import glob
import statistics
from statistics import mean
from datetime import datetime 

CWD = os.getcwd()
dirs = []           #List of all directories with (ENST, root) format, skip position 0
for x in os.walk(CWD):
    y = x[0]
    z = y.replace(CWD+"/","")    #Need to fix this replace command depending on whether its run on windows or linux
    dirs.append((z,x[0]))

enst_dic = {}
corrupted_dirs = []
for x in range(1,len(dirs)):    #parses through the list of directory paths to check if files are present and intact. It will then extract relevent data
    os.chdir(dirs[x][1])
    outfile = glob.glob("*.out")
    esfile = glob.glob("*.gff3")
    if len(outfile) == 0:
        corrupted_dirs.append(dirs[x][0])
        print(dirs[x][0])
        continue
    if len(esfile) == 0:
        corrupted_dirs.append(dirs[x][0])
        print(dirs[x][0])
        continue
    with open(outfile[0], "r") as of:   #Opens Scnfold out file and parses data
        lines = of.readlines()
        dg = []
        zs = []
        zs_one = 0
        zs_two = 0
        windows = (len(lines)-1)
        for y in range(1,len(lines)):
            line = lines[y].split()
            dg.append(float(line[3]))
            zs.append(float(line[4]))
            if float(line[4]) <= -2:
                zs_two += 1
                zs_one += 1
            elif float(line[4]) <= -1:
                zs_one += 1
            else:
                pass
        if len(dg) == 0:
            corrupted_dirs.append(dirs[x][0])
            print(dirs[x][0])
        else:
            avg_dg = mean(dg)
            avg_zs = mean(zs)
            seq_len = (windows+119)
            percent_zs1 = (100*zs_one/windows)
            percent_zs2 = (100*zs_two/windows)
    with open(esfile[0], "r") as es:    #opens the extracted structure file and pulls out motif info
        lines = es.readlines()
        motifs = len(lines)
        motif_per_len = motifs/seq_len
    if dirs[x][0] in corrupted_dirs:
        pass
    else:
        enst_dic[dirs[x][0]] = (avg_dg,avg_zs,windows,zs_one,percent_zs1,zs_two,percent_zs2,seq_len,motifs,motif_per_len)   #put info into a dictionary
os.chdir(CWD)
current_date = datetime.now().date()
current_time = datetime.now().time()
outfile1 = f"Transcriptome_metrics_{current_date}_{current_time}.txt"   #grabs date and time for output file tracking

with open(outfile1, "w") as tm:     #writes data to an output file
    tm.write("ENST\tAvg_dG\tAvg_ZS\t# of Windows\t# of ZS windows <=-1\t% of ZS windows <=-1\t# of ZS windows <=-2\t% of ZS windows <=-2\tSequence Length\t# of Motifs\tMotifs/Sequence Length\n")
    for x in range(1,len(dirs)):
        if dirs[x][0] in corrupted_dirs:
            pass
        else:
            enst = dirs[x][0]
            data = enst_dic[enst]
            tm.write(f"{enst}\t{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}\t{data[8]}\t{data[9]}\n")
outfile2 = f"Transcriptome_Metrics_corrupted_directories_{current_date}_{current_time}.txt"
with open(outfile2, "w") as cd:     #writes a file of corrupted directories which were missing needed data files
    for obj in corrupted_dirs:
        cd.write(obj+"\n")

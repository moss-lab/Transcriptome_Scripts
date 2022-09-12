#ct_sensitivity_PPV_120 
"""
By Collin A. O'Leary, Moss Lab, Iowa State University
Sensitivity and PPV comparison between two .ct files when files are run on the same input fasta, but vary in analysis parameters (e.g. inclusion of probing constrants or known reference models)
Command line to run the code: python ct_sensitivity_PPV.py ref_file.ct predicted_file.ct [output_name]
"""
import sys
import os

r_bpij = []				#A list which will hold the paired i and j values for the reference .ct file, the 'r' prefix is used throughout on variables relating to the reference file
p_bpij = []				#A list which will hold the paired i and j values for the predicted .ct file, the 'p' prefix is used throughout on variables relating to the predicted file

filename1 = sys.argv[1]				#The first system argument, this needs to be the reference file (i.e. the known structure model)
filename2 = sys.argv[2]			    #The second system argument, this is the predicted file for comparison to the reference 
#outputname = sys.argv[3]            #Specify an output file name

consistent_pairings = 0             #A running number of consistent basepairs between the files
conflicting_pairings = 0            #A running number of conflicting basepairs between the files

with open(filename1 , 'r') as rbp:			#Open reference file in read mode and refer to it as short name, then splits the file into lines
	r_lines = rbp.readlines()[1:]
    
with open(filename2 , 'r') as pbp:			#Open predicted file in read mode and refer to it as short name, then splits the file into lines
	p_lines = pbp.readlines()[1:]
 
for r_line in r_lines:                      #This loop splits the r-file lines and then iterates through it, appending the i and j nt values to a list
    data = r_line.split()
    r_bpij.append((data[0],data[4]))

for p_line in p_lines:                      #This loop splits the p file lines and then iterates through it, appending the i and j nt values to a list
    data = p_line.split()
    p_bpij.append((data[0],data[4]))
   
if len(r_bpij) >= len(p_bpij):              #This determines the shorter list and sets that length to t (which is used in range below). In some instances, ct files have an additional return at the end of the file. This accounts for that length discrepancy
    t = len(p_bpij)
else:
    t = len(r_bpij)
    
for x in range(0, t):                                                   #A loop which iterates a number of times equal to the length of the shortest list
    i = int(r_bpij[x][0])
    j = int(r_bpij[x][1])
    bp_span = abs(j-i)
    
    if bp_span <= 120:                                              #This if loop checks whether the base pair span is less than or equal to 120, if les, then it performs calculations, if over it passes
        
        if int(r_bpij[x][1]) != 0 and r_bpij[x] == p_bpij[x]:         #This checks whether a tuple holding i and j values for the reference file is paired equal to that of the predicted file
            consistent_pairings += 1                                        #If consistent, 1 is added to a rolling number which serves as a count
        elif int(r_bpij[x][1]) != 0 and r_bpij[x] != p_bpij[x]:
            conflicting_pairings += 1                                       #If conflicting, 1 is added to a rolling number which serves as a count
        else:                                                               #Unpaired nts are skipped
            pass
    else:
        
        pass
r_paired = 0                    #A variable which keeps track of the number of paired nts in the reference .ct file
r_unpaired = 0                  #A variable which keeps track of the number of unpaired nts in the reference .ct file
for obj in r_bpij:              #A loop that iterates through the reference bp list and checks whether or not the nt is paired and reports to the appropriate variable
    i = int(obj[0])
    j = int(obj[1])
    bp_span = abs(j-i)
    if bp_span <= 120:          #Only the base pairs which are 120 nts or less are addeded here as those were the only positions assessed above
        if int(obj[1]) == 0:
            r_unpaired += 1        
        else:
            r_paired += 1
    else:
        pass
p_paired = 0                    #A variable which keeps track of the number of paired nts in the predicted .ct file
p_unpaired = 0                  #A variable which keeps track of the number of unpaired nts in the predicted .ct file
for obj in p_bpij:              #A loop that iterates through the predicted bp list and checks whether or not the nt is paired and reports to the appropriate variable
    i = int(obj[0])
    j = int(obj[1])
    bp_span = abs(j-i)
    if bp_span <= 120:          #Only the base pairs which are 120 nts or less are addeded here as those were the only positions assessed above
        if int(obj[1]) == 0:
            p_unpaired += 1
        else:
            p_paired += 1

sensitivity = consistent_pairings/r_paired      #calculates the sensitivity

PPV = consistent_pairings/p_paired              #Calculates PPV


#with open(f"{outputname}.txt", "w") as file:      #Creates and opens an outputfile with specified output name
   # first_header = 'Reference File:'              #A set of different headers to be used below
    #second_header = 'Predicted File:'
    
    #Here headers and corresponding data are written to the output file
    #file.write(f'{first_header}\t{filename1}\n{second_header}\t{filename2}\nSensitivity=\t{sensitivity}\nPPV=\t{PPV}')
#file.close()
#print(consistent_pairings,r_paired,p_paired)
print(f'sensitivity=\t{sensitivity}\tPPV=\t{PPV}\tfor\t{filename1}\tand\t{filename2}')
            #Sensitivity and PPV are printed but can be written to an output file if the with open loop above is unmuted













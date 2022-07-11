# Transcriptome_Scripts

This README will detail how to use the python scripts used for analyses in "Expansion of the RNAStructuromeDB to include secondary structural data spanning the human protein coding transcriptome".
additional descriptions of the scripts can be found in the script headers.

The script HTP_dG_ZScore.pl was used to analyze the sequences of cis-regulatory elements in the Rfam database and determine their z-scores. 
To successfully run the script, it it needs to be within a directory which contains the target fasta file to be analyzed. The input file is a single fasta with any number of entries and the input value is the number of randomizations needed.
The script functions in a similar way to ScanFold-Scan to generate an output file containing a MFE, z-score, and p-value for each entry in the fasta file.

The script transcriptome_metrics was used to analyze ScanFold generated thermodynamic metrics of every transcript in the dataset.
To successfully run the script, it needs to be within a directory which contains all ScanFold output sub-directories that will be analyzed. No input files are needed.
The script will iteratively enter each sub-directory and it will find the .out file and .gff3 file produced by ScanFold.
The script will then pull relevant data out of these files to analyze and organize for output.
The script will generate two output .txt files. One will have transcriptome wide metrics formatted into a table which can be copied into an excel workbook for further analyses.
The second .txt file will contain a list of "corrupted" directories. These are directories which are missing key files or data points needed for analyses.
Ideally, the corrupted directories will be zero and this is merely a quality control step to ensure all directories were analyzed.

The script regional_zavg.py was used to analyze the per-nucleotide z-scores across a transcripts 5'UTR, CDS, and 3'UTR.
To successfully run this script, it needs to be within a directory which contains all ScanFold output sub-directories that will be analyzed. No input files are needed.
Additionally, the script needs a GFF3 file which contains regional transcript positions with their associated ENST IDs.
The script will build a dictionary using the ENST IDs and the asoociated regional position data.
Each transcript in the analyzed dataset is then compared to the dictionary and using the per-NT_average_z-score WIG file, the per-NT z-scores are binned to each region and averaged.
The script will generate two output .txt files. One will have z-score data formatted into a table which can be copied into an excel workbook for further analyses.
The second .txt file will contain a list of "corrupted" directories. These are directories which are missing key files or data points needed for analyses.
Ideally, the corrupted directories will be zero and this is merely a quality control step to ensure all directories were analyzed.

The differnetial_expression_metrics.py script was used to parse and bin out the transcriptome wide gnerated data to sub lists of specific genes.
To successfully run this script two input files are needed. The first is the output file generated from the transcriptome_metrics.py script.
The second input .txt file is a single column list of ENST ids. In our analyses, these were list of different gene expression groups.
The script will pull out the relevent data for each ENST ID from the transcriptome_metrics output and re-write the data to a new output file.

The script cm_power_parser.py was used to organize resulting out put from the cm-builder pipeline.
To successfully run this script, execute it in a directory which contains POWER files, resulting from R-scape and the cm-builder pipeline. No input files are needed.
The script will grab all available POWER files and pull out the data for base pairs which have evidence of statistically significant covariation.
The script will then sort each covarying base pair into bins (0-0.1, 0.1-0.25, >0.25 by the value of its associated power metric.
The script will produce a .txt output file with the data organized into lines for each associated power file analyzed.
The output can be viewed as a .txt file or imported into excel for more analyses.

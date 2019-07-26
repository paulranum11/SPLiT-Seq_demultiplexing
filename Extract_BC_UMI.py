#!/usr/bin/env python

import sys
import os
import argparse

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-F', '--inputFastqF', required=True, help='Input Forward Fastq File')
parser.add_argument('-R', '--inputFastqR', required=True, help='Input Reverse Fastq File')
args = parser.parse_args()

# Extract the cell barcode from the input file name
filename = str(str(args.inputFastqF))
filename_parts = filename.split(".")
Cell_Barcode = str(filename_parts[0])
Cell_Barcode_nodash_parts = Cell_Barcode.split("-")
Cell_Barcode_nodash = str(Cell_Barcode_nodash_parts[0] + Cell_Barcode_nodash_parts[1] + Cell_Barcode_nodash_parts[2])

# Create a forward output file that we will write to later
Fwd_file_name = str(Cell_Barcode + "_1.fastq")

# Make empty dictionaries
Fwd_Name_Dict = {}
Fwd_Reads_Dict = {}
Fwd_Quality_Dict = {}

Rev_Name_Dict = {}
Rev_Reads_Dict = {}
Rev_Quality_Dict = {}

Barcodes_Dict = {}
UMI_Dict = {}

# Make a position finder string. This is a static sequence in barcode 3 directly adjacent to the 8bp SPLiT-Seq barcode. 
# If the barcode sequences are changed from the standard published barcodes then this sequence needs to be modified.
PF_String = "GTGGCCG"

# Extract each line from input Forward fastq file as a string and add it to a dictionary keyed by read number.
with open(args.inputFastqF, "r") as infile:
    line_ct = 0
    for line in infile:
        if (line_ct % 4 == 0):
            Fwd_Name_Dict[line_ct]=str(line[0:].rstrip())
        if (line_ct % 4 == 1):
            Fwd_Reads_Dict[line_ct]=str(line[0:].rstrip())
        if (line_ct % 4 == 3):
            Fwd_Quality_Dict[line_ct]=str(line[0:].rstrip())
        line_ct += 1

# Extract each line from input Reverse fastq file as a string and add it to a dictionary keyed by read number.
with open(args.inputFastqR, "r") as infile:
    line_ct = 0
    for line in infile:
        if (line_ct % 4 == 0):
            Rev_Name_Dict[line_ct]=str(line[0:].rstrip())
        if (line_ct % 4 == 1):
            Rev_Reads_Dict[line_ct]=str(line[0:].rstrip())
        if (line_ct % 4 == 3):
            Rev_Quality_Dict[line_ct]=str(line[0:].rstrip())
        line_ct += 1

# Find a static consensus sequence and use it's position to find and extract UMIs. 
# If a consensus sequence is not found take the first 10bp of the read.
for key in Rev_Reads_Dict.keys():
    read = Rev_Reads_Dict[key]
    if PF_String in read[0:50]:
        strPosition1 = read.find(PF_String,0,30)
        UMI=read[int(strPosition1 - 18):int(strPosition1 - 8)]
        if (len(UMI) == 10):
            UMI_Dict[key]=UMI
        else:
            UMI_Dict[key]=read[0:10]
    else:
        UMI_Dict[key]=read[0:10]

# Write the Name, Read, and Quality scores back to the F (cDNA sequence containing) fastq file. Append the Cell and UMI barcode information to the readID.
with open(Fwd_file_name, "w+") as f:
    for key in Fwd_Reads_Dict.keys():
        ReadName = Fwd_Name_Dict[int(key - 1)]
        ReadName_parts = ReadName.split()
        print(ReadName_parts[0] + "_" + Cell_Barcode_nodash + "_" + UMI_Dict[key] + " " + ReadName_parts[1], file = f)
        print(Fwd_Reads_Dict[key], file = f)
        print("+", file = f)
        print(Fwd_Quality_Dict[int(key + 2)], file = f)

###############
# Next Steps  #
###############

#### The following bash command can be used to run this script on all the extracted fastq files.
# parallel python3 Extract_BC_UMI.py -F {} -R {}-MATEPAIR ::: $(ls *.fastq)

#### Next, a merged fastq file can be created from the reads containing the cDNA information.  
# cat *_1.fastq > MergedCells

#### At this point all the individual fastq reads can be discarded
# parallel rm {} ::: $(*fastq*)

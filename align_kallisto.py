#!/usr/bin/env python

import sys
import os
import argparse

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-F', '--inputFastqF', required=True, help='Input Fastq File')
args = parser.parse_args()

# Create empty dictionaries
Input_Dict = {}
Name_Dict = {}
Reads_Dict = {}
Quality_Dict = {}

# Name output file
Input_Dict[1] = args.inputFastqF
file = Input_Dict[1].split(".")
file = str(file[0] + ".umi")

# Define eprint function to print to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Read in input .fastq and store lines as a python dictionary
with open(args.inputFastqF, "r") as infile:
    line_ct = 0
    for line in infile:
        if (line_ct % 4 == 0):
            Name_Dict[line_ct]=str(line[1:])
        if (line_ct % 4 == 1):
            #Reads.append(str(line[1:]))
            Reads_Dict[line_ct]=str(line[1:].rstrip())
        if (line_ct % 4 == 3):
            Quality_Dict[line_ct]=str(line[1:])
        line_ct += 1

with open(file, "w+") as f:
	for key in Name_Dict.keys():
		split1 = Name_Dict[key].split("_")
		split2 = split1[1].split(" ")
		print(split2[0], end = "\n", file=f)
sys.exit()

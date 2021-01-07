#!/usr/bin/env python
# This script is intended to replace the bash wrapper script "splitseqdemultiplex_0.2.1.sh"
# It will facilitate demultiplexing with parallel processing exclusively using python

import DemultiplexUsingBarcodes_New_V2
from joblib import Parallel, delayed
import multiprocessing
import os
import Splitseq_fun_lib
import time
import sys
import argparse
import itertools
import json


#numcores = 4
#fastq_F = "fastq_F.fastq"

# Read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--numCores', required=True, help='The number of availible CPUs for parallelization.')
parser.add_argument('-e', '--errors', required=True, help='The max number of errors permissible per barcode "-e = 1" is recommended.')
parser.add_argument('-m', '--minReads', required=True, help='The minimum number of reads per cell for a cell to be retained.')
parser.add_argument('-1', '--round1Barcodes', required=True, help='output .fastq')
parser.add_argument('-2', '--round2Barcodes', required=True, help='output .fastq')
parser.add_argument('-3', '--round3Barcodes', required=True, help='output .fastq')
parser.add_argument('-f', '--fastqF', required=True, help='filepath to the forward .fastq file')
parser.add_argument('-r', '--fastqR', required=True, help='filepath to the reverse .fastq file')
parser.add_argument('-o', '--outputDir', required=True, help='name of the output directory')
parser.add_argument('-t', '--targetMemory', required=True, help='filepath to the reverse .fastq file')
parser.add_argument('-a', '--align', required=True, help='filepath to the reverse .fastq file')
parser.add_argument('-x', '--starGenome', required=True, help='filepath to the reverse .fastq file')
parser.add_argument('-y', '--starGTF', required=False, help='filepath to the reverse .fastq file')
parser.add_argument('-s', '--geneAnnotationSAF', required=False, help='filepath to the reverse .fastq file')
parser.add_argument('-b', '--numReadsBin', required=True, help='the number of reads to be processed before results are flushed to disc')
parser.add_argument('-p', '--performanceMetrics', required=True, help='Enter True or False to turn on or off reporting of the number of demultiplexed cells') 
parser.add_argument('-l', '--lengthFastq', required=True, help='the length (number of lines) inthe provided input fastqR file. This can be obtained using the "wc -l fastqR" command on linux systems.')
args = parser.parse_args()

#####################################
# STEP1: Demultiplex by position    #
#####################################
print("Starting Step1: Splitting input fastq files.")
#starttime = time.time()

# Run the split_fastq_fun function to split the input fastq into bins for parallel processing.
Splitseq_fun_lib.split_fastqF_fun(int(args.numCores), args.fastqF, int(args.lengthFastq))
Splitseq_fun_lib.split_fastqR_fun(int(args.numCores), args.fastqR, int(args.lengthFastq))

# Execute demultiplexing script using parallel processing
Parallel(n_jobs = int(args.numCores))(delayed(DemultiplexUsingBarcodes_New_V2.demultiplex_fun)(inputFastqF = str("split_fastq_F_" + str(i)), 
    inputFastqR = str("split_fastq_R_" + str(i)),
    outputDir = args.outputDir,
    numReadsBin = int(args.numReadsBin),
    errorThreshold = int(args.errors),
    performanceMetrics = True,
    readsPerCellThreshold = 10,
    verbose = True) for i in range(int(args.numCores)))

###########################################################
# STEP2: Extract Performance Metrics and Apply Filtering  #
###########################################################
# Calculate number of barcodes metrics from the generated .fastq files
DemultiplexUsingBarcodes_New_V2.calc_demux_results(outputDir = args.outputDir, performanceMetrics = bool(args.performanceMetrics), readsPerCellThreshold = args.minReads)




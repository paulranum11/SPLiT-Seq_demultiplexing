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
parser.add_argument('-t', '--numCores', required=True, help='The number of availible threads for parallelization.')
parser.add_argument('-e', '--errors', required=True, help='The max number of errors permissible per barcode "-e = 1" is recommended.')
parser.add_argument('-m', '--minReads', required=True, help='The minimum number of reads per cell for a cell to be retained.')
parser.add_argument('-1', '--round1Barcodes', required=True, help='output .fastq')
parser.add_argument('-2', '--round2Barcodes', required=True, help='output .fastq')
parser.add_argument('-3', '--round3Barcodes', required=True, help='output .fastq')
parser.add_argument('-f', '--fastqF', required=True, help='filepath to the forward .fastq file')
parser.add_argument('-r', '--fastqR', required=True, help='filepath to the reverse .fastq file')
parser.add_argument('-o', '--outputDir', required=True, help='name of the output directory')
parser.add_argument('-a', '--align', required=True, help='filepath to the reverse .fastq file')
parser.add_argument('-x', '--starGenome', required=False, help='filepath to the reverse .fastq file')
parser.add_argument('-s', '--geneAnnotationSAF', required=False, help='filepath to the reverse .fastq file')
parser.add_argument('-b', '--numReadsBin', required=True, help='the number of reads to be processed before results are flushed to disc')
parser.add_argument('-p', '--performanceMetrics', required=True, help='Boolean (True or False) to turn on or off reporting of the number of demultiplexed cells') 
parser.add_argument('-l', '--lengthFastq', required=True, help='the length (number of lines) inthe provided input fastqR file. This can be obtained using the "wc -l fastqR" command on linux systems.')
parser.add_argument('-n', '--nanopore', required=False, help='Boolean (True or False) to control use of minimap2 as the aligner to facilitate alignment of nanopore long read data.')
parser.add_argument('-g', '--minimap2Genome', required=False, help='Enter filepath to minimap2 genome if the -n True flag is used')
args = parser.parse_args()

#####################################
# STEP1: Demultiplex by position    #
#####################################
print("Starting Step1: Splitting input fastq files.")
#starttime = time.time()

# Create the output directory
#if not os.path.exists(args.outputDir):
#    os.makedirs(args.outputDir)

# Run the split_fastq_fun function to split the input fastq into bins for parallel processing.
#Splitseq_fun_lib.split_fastqF_fun(int(args.numCores), args.fastqF, int(args.lengthFastq))
#Splitseq_fun_lib.split_fastqR_fun(int(args.numCores), args.fastqR, int(args.lengthFastq))

# Execute demultiplexing script using parallel processing
#Parallel(n_jobs = int(args.numCores))(delayed(DemultiplexUsingBarcodes_New_V2.demultiplex_fun)(inputFastqF = str("split_fastq_F_" + str(i)), 
#    inputFastqR = str("split_fastq_R_" + str(i)),
#    outputDir = args.outputDir,
#    numReadsBin = int(args.numReadsBin),
#    errorThreshold = int(args.errors),
#    performanceMetrics = bool(args.performanceMetrics),
#    readsPerCellThreshold = 10,
#    verbose = True) for i in range(int(args.numCores)))

###########################################################
# STEP2: Extract Performance Metrics and Apply Filtering  #
###########################################################
# Calculate number of barcodes metrics from the generated .fastq files
#DemultiplexUsingBarcodes_New_V2.calc_demux_results(outputDir = args.outputDir, performanceMetrics = bool(args.performanceMetrics), readsPerCellThreshold = args.minReads)


###########################################################
# STEP3: Remove Intermediate Files                        #
###########################################################
#for i in range(int(args.numCores)):
#    Splitseq_fun_lib.remove_file_fun(str("./split_fastq_F_" + str(i)))
#    Splitseq_fun_lib.remove_file_fun(str("./split_fastq_R_" + str(i)))

#Splitseq_fun_lib.remove_file_fun("./position_learner_fastqr.fastq")
#Splitseq_fun_lib.remove_file_fun(str("./" + args.outputDir + "/MergedCells_1.fastq"))
#Splitseq_fun_lib.remove_dir_fun("./__pycache__")


##########################################################
# STEP4: Perform Mapping                                 # 
##########################################################
counter = 0
while counter <= 1:
    if bool(args.align) == True:
        if bool(args.nanopore) == True:
            Splitseq_fun_lib.run_minimap2_alignment_fun(numCores = args.numCores, resultsDir = args.outputDir, minimap2Genome = args.minimap2Genome)
            Splitseq_fun_lib.run_featureCounts_SAF_fun(numCores = args.numCores, countsFile = args.geneAnnotationSAF, resultsDir = args.outputDir)
            Splitseq_fun_lib.run_samtools_fun(resultsDir = args.outputDir)
            Splitseq_fun_lib.run_umi_tools_fun(resultsDir = args.outputDir)
            counter += 1
            break
        else: 
            continue
        
        print("Alignment setting = True")
        Splitseq_fun_lib.run_star_alignment_fun(numCores = args.numCores, starGenome = args.starGenome, resultsDir = args.outputDir)
        Splitseq_fun_lib.run_featureCounts_SAF_fun(numCores = args.numCores, countsFile = args.geneAnnotationSAF, resultsDir = args.outputDir)
        Splitseq_fun_lib.run_samtools_fun(resultsDir = args.outputDir)
        Splitseq_fun_lib.run_umi_tools_fun(resultsDir = args.outputDir)
        counter += 1
    else:
        print("Alignment setting = False, exiting run")

##########################################################
# STEP5: Final Clean Up
##########################################################
Splitseq_fun_lib.remove_dir_fun("__pycache__")

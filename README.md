# Latest Updates
11-25-2018 - Speed was dramatically improved through modifications to the matepair identification step.

# Splitseq_demultiplexing
An unofficial demultiplexing strategy for SPLiT-seq RNA-Seq data.  This approach DOES NOT conform to the exact specifications reported in the SPLiT-Seq paper. It will produce one .fastq file per individual cell sample as defined by their unique barcode configuration.  

# System requirements
This script has been tested on a linux cluster running Linux 3.10.0-514.2.2.e17.x86_64 

This script is written in bash and should be portable across a variety of linux systems running the bash shell.

In order to run this software you must install the following dependency packages.

PYTHON3 needs to be set as your default python installation

GNU parallel: https://www.gnu.org/software/parallel/

UMI tools: https://github.com/CGATOxford/UMI-tools

# Getting Started
Download this git repository .zip file or clone this repository using `git clone`. The downloaded directory will contain three (Round1, Round2, and Round3) barcode files as well as a small example dataset derrived from the 100_CNS_nuclei dataset GEO assession: GSM3017260 (SRR6750041).  The full sized datasets can be downloaded from the following European Nucleotide Archive address https://www.ebi.ac.uk/ena/data/view/SRR6750041

The executable file is called `splitseqdemultiplex.sh` it is written in for bash and can be called using `bash splitseqdemultiplex.sh (options)`

# Options
-n | --numcores # specifies the number of cores you would like to use to parallelize your run.

-m | --minreads # specifies the minimum number of reads required for a cell to be retained. The default is set to 10.

-1 | --round1barcodes # specifies name of the file containing the barcodes you would like to use for round1. These should be provided as a separate file.  See the provided example for formatting reference.

-2 | --round2barcodes # specifies name of the file containing the barcodes you would like to use for round2. These should be provided as a separate file.  See the provided example for formatting reference.

-3 | --round3barcodes # specifies name of the file containing the barcodes you would like to use for round3. These should be provided as a separate file.  See the provided example for formatting reference.

-f | --fastqF # filepath to the Forward input .fastq file. 

-r | --fastqR # filepath to the Reverse input .fastq file.

-o | --outputdir # filepath to the desired output directory.

Users may increase the speed of the run by allocating additonal cores using -n and increasing the minimum number of reads required for each cell using -m.  Default values for -1 -2 and -3 are the barcodes provided in the splitseq_demultiplexing download: `Round1_barcodes_new3.txt`, `Round2_barcodes_new3.txt` and `Round3_barcodes_new3.txt`.  Default values for `-f` and `-r` are the provided example .fastq files.  The default output directory is `results`

# Example
The following is an example command that will run splitseqdemultiplex.sh using the provided example datasets.

`bash splitseqdemultiplex.sh -n 4 -m 10 -1 Round1_barcodes_new3.txt -2 Round2_barcodes_new3.txt -3 Round3_barcodes_new3.txt -f SRR6750041_1_smalltest.fastq -r SRR6750041_2_smalltest.fastq -o results`

# Benchmarking
Updated: 11-26-2018

Benchmarking was performed on a large previously published ~16Gb .fastq dataset found here https://www.ebi.ac.uk/ena/data/view/SRR6750041. `splitseqdemultiplex.sh` was run on four cores of a linux system. 

STEP1 (Demultiplexing): Time elapsed = 4hrs 35min 20sec

STEP2 (Matepair Finding): Time elapsed = 

STEP3 (UMI Extraction): Time elapsed =

Total: Time elapsed = 

# Notes and caution
This tool is under development. No warranty is implied and accurate function is NOT guarenteed.  

If you would like to contribute to this tool please help us make it better! 

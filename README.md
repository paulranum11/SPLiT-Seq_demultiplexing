# Splitseq_demultiplexing
An unofficial demultiplexing strategy for SPLiT-seq RNA-Seq data.  This approach DOES NOT conform to the exact specifications reported in the SPLiT-Seq paper. It will produce one .fastq file per individual cell sample as defined by their unique barcode configuration.  

# System requirements
This script has been tested on a linux cluster running Linux 3.10.0-514.2.2.e17.x86_64 

This script is written in bash and should be portable across a variety of linux systems running the bash shell.

In order to run this software you must install the following dependency packages.

GNU parallel: https://www.gnu.org/software/parallel/

UMI tools: https://github.com/CGATOxford/UMI-tools

# Usage
Download this git repository .zip file or clone this repository using `git clone`. The downloaded directory will contain three (Round1, Round2, and Round3) barcode files as well as a small and medium sized example dataset derrived from the 100_CNS_nuclei dataset GEO assession: GSM3017260 (SRR6750041).  The full sized datasets can be downloaded from the following European Nucleotide Archive address https://www.ebi.ac.uk/ena/data/view/SRR6750041

To run the script type `bash splitseqdemultiplex.sh`

Users may increase the speed thorugh parallelization of the matepair and UMI finding parts of the script by increasing the number of cores availible.  This can be set by entering the number of cores you wish to use in the `numcores="4"` parameter of the splitseqdemultiplex.sh file. Simply replace 4 with your desired number of cores.  If you are working on a cluster you can submit this script to your queue by adding a standard header to the file and submitting using `qsub`.

# Notes and caution
This tool is actively in development no warranty is implied and accurate function is NOT guarenteed.  

If you would like to contribute to this tool please help us make it better! 

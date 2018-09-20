# splitseq_demultiplexing
An unofficial demultiplexing strategy for SPLiT-seq RNA-Seq data.  This approach DOES NOT conform to the exact specifications reported in the SPLiT-Seq paper. It will produce as output one .fastq file per individual cell sample as defined by their barcodes.  

# system requirements
This script has been tested on a linux cluster running Linux 3.10.0-514.2.2.e17.x86_64 

This script is written in bash and should be portable across a variety of linux systems running the bash shell.

In order to run this software you must install the following dependency packages.

GNU parallel: https://www.gnu.org/software/parallel/

UMI tools: https://github.com/CGATOxford/UMI-tools

# Usage
To use this script download this git repository .zip file or clone this repository using `git clone`.



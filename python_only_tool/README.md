README.md

This directory contains a development version of an All-python splitseq_demultiplexing pipeline.

###Dependencies: 
1. Python3

2. Python3 packages: 
	- joblib
	- multiprocessing
	- os
	- sys
	- argparse
	- itertools

3. For alignment and genes per cell counts table generation you need a linux based system and the following software 
	- STAR (and a STAR genome index) `https://github.com/alexdobin/STAR`
	- featureCounts `http://subread.sourceforge.net/`	
	- umi-tools `https://github.com/CGATOxford/UMI-tools`
	- samtools `https://github.com/samtools/samtools`

###Arguments: \
	- '-n', '--numCores', required=True, help=The number of availible threads for parallelization. \
	- '-e', '--errors', required=True, help=The max number of errors permissible per barcode "-e = 1" is recommended. \
	- '-m', '--minReads', required=True, help=The minimum number of reads per cell for a cell to be retained. \
	- '-1', '--round1Barcodes', required=True, help=The path to the provided round1Barcodes text file. \
	- '-2', '--round2Barcodes', required=True, help=The path to the provided round2Barcodes text file. \
	- '-3', '--round3Barcodes', required=True, help=The path to the provided round3Barcodes text file. \
	- '-f', '--fastqF', required=True, help=The path to the forward .fastq file \
	- '-r', '--fastqR', required=True, help=The path to the reverse .fastq file \
	- '-o', '--outputDir', required=True, help=The name of the output directory (this directory will be created by the script). \
	- '-a', '--align', required=True, help=Boolean (True or False), indicating if you want to perform STAR alignment and genes per cell counts matrix  generation. \
	- '-x', '--starGenome', required=True, help=The path to the STAR genome. \
	- '-s', '--geneAnnotationSAF', required=False, help=The path to a SAF annotation file corresponding to your STAR genome. \
	- '-b', '--numReadsBin', required=True, help=the number of reads to be processed before results are written to disc (if you are unsure try 100000). \
	- '-p', '--performanceMetrics', required=True, help=boolean (True or False) to turn on or off reporting of the number of demultiplexed cells. This step is required to filter the results .fastq file based on the `--minReads` threshold provided above. \
	- '-l', '--lengthFastq', required=True, help=the length (in number of lines) of the input fastqR file. This can be obtained using the "wc -l fastqR" command on linux systems. \
        - '-k', '--positionDetection', --positionDetection', required=False, dest='positionDetection', action='store_false', help=provide only the -k flag to turn off automatic barcode position detection and instead use the default barcode positions. 

###Example Use: 
python splitseqdemultiplex_0.2.2.py \
	-n 4 # For a full scale run on large input files we recommend that you use as many availible cores as possible to increase run speed. \
	-e 1 \
	-m 2 # For a full scale run we recommend that you use a minimum reads per cell threshold of 200 or greater to reduce need for downstream filtering and to help constrain the size of the results counts.tsv.gz file \
	-1 Round1_barcodes_new5.txt \
	-2 Round2_barcodes_new4.txt \
	-3 Round3_barcodes_new4.txt \
	-f SRR6750041_1_smalltest.fastq # For a full scale run please provide the filepaths to your forward (_1_) and reverse (_2_) fastq file. \
	-r SRR6750041_2_smalltest.fastq \ 
	-o output \
	-a False \
	-x /path/to/my/star/genome \
	-b 100000 \
	-p True \
	-l 3000000	


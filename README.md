# SPLiT-Seq_demultiplexing_0.1.1
This tool was created to provide an open source, portable solution for demultiplexing SPLiT-Seq RNA-Seq datasets. SPLiT-Seq_demultiplexing has two core versions:
1. `--version merge` which produces one .fastq file in which CellID and UMI information is appended to the readID. 
2. `--version split` which produces one .fastq file for each single-cell identified. Output .fastq files are named using the identified barcode combination and UMIs are appended to the readID.

During demultiplexing each "cell" is defined by its unique configuration of SPLiT-Seq round1-3 barcodes. 

Optional alignment `--align`, gene assignment, and counts per gene (per cell) table generation functionality is included.  


# System Requirements
This script has been tested on a linux cluster running Linux CentOS 3.10.0-514.2.2.e17.x86_64 and on a MacBook Pro running macOS High Sierra v10.13.6. NOTE: on macOS systems the options do not work. This problem can be resolved by hardcoding the options you want inside the splitseqdemultiplexing_0.1.0.sh under `set default inputs` or by installing GNU getopt. 

This script is written in bash and python3 and should be portable across a variety of linux systems running the bash shell.

In order to run this software you must install the following dependency packages.

- Python3 needs to be installed on your system. Often the executable name of python3 can vary... for example it may appear as `python` or as `python3`. This script requires that the executable be `python`. 
- Python3 packages: `math`, `os`, `psutil`, `argparse`, `sys`, `datetime`, `itertools`, `re`
- GNU parallel: `https://www.gnu.org/software/parallel/`
- UMI-tools: `https://github.com/CGATOxford/UMI-tools`
- STAR: `https://github.com/alexdobin/STAR` in order to run the optional aligment step the STAR rna-seq aligner is required along with an appropriate STAR genome index.
- featureCounts: `http://subread.sourceforge.net/` 
- Samtools: `https://github.com/samtools/samtools` 


# Getting Started
Download this git repository .zip file or clone this repository using `git clone`. The downloaded directory will contain three (Round1, Round2, and Round3) barcode files as well as a small example dataset derrived from the 100_CNS_nuclei dataset GEO accession: GSM3017260 (SRR6750041).  The full sized datasets can be downloaded from the following European Nucleotide Archive address https://www.ebi.ac.uk/ena/data/view/SRR6750041

The executable file is called `splitseqdemultiplex_0.1.0.sh` it is written in bash and can be called using `bash splitseqdemultiplex_0.1.0.sh (options)`


# Options
-n | --numcores # specifies the number of cores you would like to use to parallelize your run.

-v | --version # specifies the version of the demultiplexing utility you would like to run. Input `merged` to output a single .fastq file with CellIDs and UMIs annotated on each read ID line of the .fastq file.  Input `split` to output a single .fastq file for every single cell sample identified in the input .fastq file with UMIs appended to the read ID line of the .fastq file.

-e | --errors # specifies the number of errors acceptable at each barcode position. The default is set to 1.

-m | --minreads # specifies the minimum number of reads required for a cell to be retained. The default is set to 10.

-1 | --round1barcodes # specifies name of the file containing the barcodes you would like to use for round1. These should be provided as a separate file.  See the provided example for formatting reference.

-2 | --round2barcodes # specifies name of the file containing the barcodes you would like to use for round2. These should be provided as a separate file.  See the provided example for formatting reference.

-3 | --round3barcodes # specifies name of the file containing the barcodes you would like to use for round3. These should be provided as a separate file.  See the provided example for formatting reference.

-f | --fastqF # filepath to the Forward input .fastq file. 

-r | --fastqR # filepath to the Reverse input .fastq file.

-o | --outputdir # filepath to the desired output directory.

-t | --targetMemory # define the memory maximum. Processed reads will be saved to memory until this memory maximum is reached.  A higher value increases the speed of the script but uses more system memory. Our recommended value is `8000` which equates to 8gb. Higher or lower values will work fine but we suggest using more if your system can support it. 

-g | --granularity # the granularity with which you want to get progress updates. Default value is `100000`.

-c | --collapseRandomHexamers # when `true` this option will collapse unique barcode combinations primed with Random Hexamers and OligoDT primers.  Because SPLiT-Seq uses both Random Hexamers and OligoDT primers with different barcodes in the same well of the ROUND1 RT step this option is set to `true` by default.

-a | --align # This is an optional argument. If not included the script will terminate after producing the output .fastq file or files. If `star` or `kallisto` are entered as inputs an alignment will be initiated. The `star` aligner must be used with the `-v merged` argument and the `kallisto` aligner may only be used with the `-v split` option. After alignment, per cell, gene level expression abundance will be calculated and a counts matrix will be produced.

-x | --starGenome # provide the path to the STAR genome index file. Note: STAR options are only relevant if you are doing STAR alignment.
 
-y | --starGTF # The format of this argument is "GTF /path/to/my/file.gtf" IMPORTANT: The letters GTF must be included before the filepath and the entire argument should be wrapped in quotes as shown in the example. The .gtf file used should be the genome annotation .gtf file that corresponds to the STAR index genome. This is not a file generated by STAR but just a GTF format file downloaded from the same place you downloaded your raw genome. IMPORTANT: Use either a .gtf index using -y or a .saf index using -s. Do not use both.

-s | --geneAnnotationSAF # The format of this argument is "SAF /path/to/my/file.gtf" IMPORTANT: The letters SAF must be included before the filepathand the entire argument should be wrapped in quotes as shown in the example. The .saf file must correspond to the genome used to construct your STAR index. SAF format annotation files are preferred for nuclei sequencing when the user wants to assign both intronic and exonic reads to genes. Prebuilt indexes are provide for Human (GRCh38) and Mouse (GRCm38) from Ensembl.org. IMPORTANT: Use either a .gtf index using -y or a .saf index using -s. Do not use both.

-f | --kallistoIndexIDX # provide the path to your kallisto index file (.idx).  Note: Kallisto options are only relevant if you are doing kallisto alignment.

-i | --kallistoIndexFASTA # provide the path to the .fasta file corresponding to your kallisto index file. 

Notes:
Users may increase the speed of the run by allocating additonal cores using -n and increasing the minimum number of reads required for each cell using -m.  Default values for -1 -2 and -3 are the barcodes provided in the splitseq_demultiplexing download: `Round1_barcodes_new3.txt`, `Round2_barcodes_new3.txt` and `Round3_barcodes_new3.txt`.  Default values for `-f` and `-r` are the provided example .fastq files.  The default output directory is `results`

# Example
The following is an example command that will run splitseqdemultiplex.sh using the provided example datasets.

`bash splitseqdemultiplex.sh -n 4 -e 1 -m 10 -1 Round1_barcodes_new5.txt -2 Round2_barcodes_new4.txt -3 Round3_barcodes_new4.txt -f SRR6750041_1_smalltest.fastq -r SRR6750041_2_smalltest.fastq -o results -t 8000 -g 100000 -c true -a star -x ~/my/path/to/starIndexDirectory/GRCm38/ -s "SAF ~/path/to/GRCm38.saf"`


# Benchmarking
Updated: Jan_16_2019

Benchmarking was performed on a previously published, ~17Gb (77,621,181 read) fastq dataset found here https://www.ebi.ac.uk/ena/data/view/SRR6750041. `splitseqdemultiplex.sh` was run on four cores of a linux (CentOS) system using `-t 8000`. A maximum of one error was permitted at each barcode position and cells containing fewer than 10 reads were discarded.

STEP1 (Demultiplexing): Time elapsed = 3hrs 49min 09sec

STEP2 (Matepair Finding): Time elapsed = 8min 35sec 

STEP3 (UMI Extraction): Time elapsed = 4hrs 52min 14sec 

NOTE: Speed is dependant on the size of the input files, the amount of memory allocated using `-t`, and the number of cores used.

# Benchmarking Output
Updated: Feb_11_2019

- 12,768 .fastq files (cells) were generated as output, containing a total of 56,256,959 reads.
- Each result contained >= 10 reads, the default minimum read cutoff. 
- The largest result .fastq file "AGCATTCGCAACCACA-ACACAGAAATCCA-TGGAACAAGTGGCC.fastq" contained 4,277,136 individual reads.
- The mean number of reads per cell was 4406.
- 34,162,777 reads were primed with OligoDT primers and 22,094,182 were primed with Random Hexamer primers.

# Latest Updates
- Jul-01-2019 - Single output .fastq format was adopted with CellID and UMI information added to readIDs. STAR alignment was implemented and UMI_tools based gene counts matrix generation was added.
- Mar-14-2019 - Support for Kallisto pseudoalignment and expression quantification was added.
- Feb-09-2019 - Support for random hexamer primers was added. When `-c` is `true` random hexamer reads will be detected and added to the cell from which they originate.  
- Jan-16-2019 - HUGE update to dramatically increase speed. STEP1 and STEP2 were completely rewritten to make use of hashing and python dictionaries. Big thanks to Charlie Whitmore for making this possible!
- Dec-18-2018 - Added support for reads containing sequencing errors. The number of permissible errors is defined by the user using -e 'number' (default = 1).
- Nov-25-2018 - Speed was dramatically improved through modifications to the matepair identification step.

# Notes and Caution
This tool is under development. No warranty is implied and accurate function is NOT guarenteed. This approach does not conform to the exact specifications reported in the SPLiT-Seq paper.

# Contributors
We welcome contributions that make this tool better! If you think you found an error or solved an issue please let us know!
Big thanks to developers who have made improvements to this tool!
- Charlie Whitmore
- Cody Raspen

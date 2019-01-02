#!/bin/bash

###############
# Example Use #
###############

#bash splitseqdemultiplex.sh \
# -n 12 \
# -e 1 \
# -m 10 \
# -1 Round1_barcodes_new3.txt \
# -2 Round2_barcodes_new3.txt \
# -3 Round3_barcodes_new3.txt \
# -f SRR6750041_1_smalltest.fastq \
# -r SRR6750041_2_smalltest.fastq \
# -o results


################
# Dependencies #
################
# Python3 must be installed and accessible as "python" from your system's path
type python &>/dev/null || { echo "ERROR python3 is not installed or is not accessible from the PATH as python"; exit 1; }

# agrep must be installed and accessible as "agrep" from your system's path
type agrep &>/dev/null || { echo "ERROR agrep is not installed or is not accessible from the PATH as agrep"; exit 1; }

# UMI_Tools must be installed and accessible from the PATH as "umi_tools"
type umi_tools &>/dev/null || { echo "ERROR umi_tools is not installed or is not accessible from the PATH as umi_tools"; exit 1; }

# parallel must be installed and accessible from the path as "parallel"
type parallel &>/dev/null || { echo "ERROR parallel is not installed or is not accessible from the PATH as parallel"; exit 1; }


###########################
### Set Default Inputs  ###
###########################

NUMCORES="4"
ERRORS="1"
MINREADS="10"
ROUND1="Round1_barcodes_new3.txt"
ROUND2="Round2_barcodes_new3.txt"
ROUND3="Round3_barcodes_new3.txt"
FASTQ_F="SRR6750041_1_smalltest.fastq"
FASTQ_R="SRR6750041_2_smalltest.fastq"
OUTPUT_DIR="results"



################################
### User Inputs Using Getopt ###
################################
# read the options
TEMP=`getopt -o n:m:1:2:3:f:r:o: --long numcores:,errors:,minreads:,round1barcodes:,round2barcodes:,round3barcodes:,fastqF:,fastqR:,outputdir: -n 'test.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -n|--numcores)
            case "$2" in
                "") shift 2 ;;
                *) NUMCORES=$2 ; shift 2 ;;
            esac ;;
        -e|--errors)
            case "$2" in
                "") shift 2 ;;
                *) ERRORS=$2 ; shift 2 ;;
            esac ;;
        -m|--minreads)
            case "$2" in
                "") shift 2 ;;
                *) MINREADS=$2 ; shift 2 ;;
            esac ;;
        -1|--round1barcodes)
            case "$2" in
                "") shift 2 ;;
                *) ROUND1=$2 ; shift 2 ;;
            esac ;;
        -2|--round2barcodes)
            case "$2" in
                "") shift 2 ;;
                *) ROUND2=$2 ; shift 2 ;;
            esac ;;
        -3|--round3barcodes)
            case "$2" in
                "") shift 2 ;;
                *) ROUND3=$2 ; shift 2 ;;
            esac ;;
        -f|--fastqF)
            case "$2" in
                "") shift 2 ;;
                *) FASTQ_F=$2 ; shift 2 ;;
            esac ;;
        -r|--fastqR)
            case "$2" in
                "") shift 2 ;;
                *) FASTQ_R=$2 ; shift 2 ;;
            esac ;;
        -o|--outputdir)
            case "$2" in
                "") shift 2 ;;
                *) OUTPUT_DIR=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done


###############################
### Write Input Args To Log ###
###############################

# Print the arguments provided as input to splitseqdemultiplex.sh
echo "splitseqdemultiplex.sh has been run with the following input arguments"
echo "numcores = $NUMCORES" 
echo "errors = $ERRORS" 
echo "minreadspercell = $MINREADS" 
echo "round1_barcodes = $ROUND1" 
echo "round2_barcodes = $ROUND2" 
echo "round3_barcodes = $ROUND3" 
echo "fastq_f = $FASTQ_F" 
echo "fastq_r = $FASTQ_R" 

# calculate the min number of lines per cell 
minlinesperfastq=$(($MINREADS * 4))
echo "minimum reads per cell set to $MINREADS" 
echo "minimum lines per cell is set to $minlinesperfastq" 


#######################################
# STEP 1: Demultiplex Using Barcodes  #
#######################################

# Add the barcode sequences to a bash array.
declare -a ROUND1_BARCODES=( $(cut -b 1- $ROUND1) )
declare -a ROUND2_BARCODES=( $(cut -b 1- $ROUND2) )
declare -a ROUND3_BARCODES=( $(cut -b 1- $ROUND3) )

# Log current time
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Current time : $now"  

# Make folder for $OUTPUT_DIR files
rm -r $OUTPUT_DIR
mkdir $OUTPUT_DIR
touch $OUTPUT_DIR/emptyfile.txt

# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP1: Demultiplex using barcodes. Current time : $now" 

# Create search function that uses awk and agrep to find approximate matches in the fastq reads.
srch() {
   awk -F ': ' 'NR==FNR {
      a[$1] = 1
      next
   }
   a[FNR] {
      print p
      print
      for (i=0; i<2 && getline > 0; i++)
         print
   }
   {
      p=$0
   }' <(agrep -"$ERRORS" -n "$2" "$1") "$1"
}
export -f srch

# Begin the set of nested loops that searches for every possible barcode. We begin by looking for ROUND1 barcodes 
for barcode1 in "${ROUND1_BARCODES[@]}";
    do
    srch $FASTQ_R "$barcode1" > ROUND1_MATCH.fastq
    find $OUTPUT_DIR/ -size 0 -delete 
    
        if [ -s ROUND1_MATCH.fastq ]
        then
            
            # Now we will look for the presence of ROUND2 barcodes in our reads containing barcodes from the previous step
            for barcode2 in "${ROUND2_BARCODES[@]}";
            do
            srch ROUND1_MATCH.fastq "$barcode2" > ROUND2_MATCH.fastq
                
                # Now we will look for the presence of ROUND3 barcodes in our reads containing barcodes from the previous step. This step is parallelized using GNU parallel
                if [ -s ROUND2_MATCH.fastq ]
                then
                    parallel -j $NUMCORES "srch ROUND2_MATCH.fastq {} > ./$OUTPUT_DIR/$barcode1-$barcode2-{}.fastq" ::: "${ROUND3_BARCODES[@]}"

                fi
            done
        fi
    done

# find and remove all files with 0 file size
find ./$OUTPUT_DIR -size 0 -print0 |xargs -0 rm --

# calculate the number of cells (.fastq files) with >= 1 read 
numfastqbeforeremoval=$(ls -tslh ./$OUTPUT_DIR | wc -l) 

# Create a function to remove .fastq files containing fewer than a user defined minimum number of reads
# Reads are multiplied by four to get the number of lines
removebylinesfunction () {
find "$1" -type f |
while read f; do
        i=0
        while read line; do
                i=$((i+1))
                [ $i -eq $minlinesperfastq ] && continue 2
        done < "$f"
        printf %s\\n "$f"
done |
xargs rm -f
}
export -f removebylinesfunction 

# Run the function to remove .fastq files containing fewer than the minimum number of lines
removebylinesfunction ./$OUTPUT_DIR

numfastqafterremoval=$(ls -tslh ./$OUTPUT_DIR | wc -l)

echo "$numfastqbeforeremoval cells were identified containing >= 1 read" 
echo "$numfastqafterremoval cells were identified containing >= $MINREADS reads, the minimum number of reads defined by the user." 

# Remove remaining round1 and round2 intermediate .fastq files
rm ROUND*


##########################################################
# STEP 2: For every cell find matching paired end reads  #
##########################################################
# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP2: Finding read mate pairs. Current time : $now" 

# Now we need to collect the other read pair. To do this we can collect read IDs from the $OUTPUT_DIR files we generated in step one.
# Generate an array of cell filenames
python3 matepair_finding.py --input $OUTPUT_DIR --fastqf $FASTQ_F --output $OUTPUT_DIR


########################
# STEP 3: Extract UMIs #
########################
# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP3: Extracting UMIs. Current time : $now" 

rm -r $OUTPUT_DIR-UMI
mkdir $OUTPUT_DIR-UMI

# Parallelize UMI extraction
{
#parallel -j $NUMCORES 'fastp -i {} -o results_UMI/{/}.read2.fastq -U --umi_loc=read1 --umi_len=10' ::: results/*.fastq
parallel -j $NUMCORES -k "umi_tools extract -I {} --read2-in={}-MATEPAIR --bc-pattern=NNNNNNNNNN --log=processed.log --read2-out=$OUTPUT_DIR-UMI/{/}" ::: $OUTPUT_DIR/*.fastq
#parallel -j $NUMCORES 'mv {} $OUTPUT_DIR_UMI/cell_{#}.fastq' ::: $OUTPUT_DIR_UMI/*.fastq
} &> /dev/null

#################################
# STEP 4: Collect Summary Stats #
#################################
# Print the number of lines and barcode ID for each cell to a file
echo "$(wc -l results-UMI/*.fastq)" | sed '$d' | sed 's/results-UMI\///g' > linespercell.txt
Rscript generate_reads_violin.r


###########################
# STEP 4: Perform Mapping #
###########################
#parallel -j $NUMCORES "STAR --runThreadN 1 
#--genomeDir /mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/GRCh38 
#--readFilesIn {} 
#--outFilterMismatchNoverLmax 0.05 
#--alignIntronMax 20000 
#--outSAMstrandField intronMotif 
#--quantMode TranscriptomeSAM GeneCounts 
#--outSAMtype BAM SortedByCoordinate" ::: $OUTPUT_DIR-UMI/*.fastq

#All finished
number_of_cells=$(ls -1 "$OUTPUT_DIR-UMI" | wc -l)
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "a total of $number_of_cells cells were demultiplexed from the input .fastq" 
echo "Current time : $now" 
echo "all finished goodbye" 

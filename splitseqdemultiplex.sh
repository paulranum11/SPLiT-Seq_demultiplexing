#!/bin/bash

#alias python='python3'

###############
# Example Use #
###############

#bash splitseqdemultiplex.sh \
# -n 12 \
# -e 1 \
# -m 10 \
# -1 Round1_barcodes_new4.txt \
# -2 Round2_barcodes_new4.txt \
# -3 Round3_barcodes_new4.txt \
# -f SRR6750041_1_smalltest.fastq \
# -r SRR6750041_2_smalltest.fastq \
# -o results \
# -t 8000 \
# -g 100000 


################
# Dependencies #
################
# Python3 must be installed and accessible as "python" from your system's path
type python &>/dev/null || { echo "ERROR python3 is not installed or is not accessible from the PATH as python"; exit 1; }

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
ROUND1="Round1_barcodes_new5.txt"
ROUND2="Round2_barcodes_new4.txt"
ROUND3="Round3_barcodes_new4.txt"
FASTQ_F="SRR6750041_1_smalltest.fastq"
FASTQ_R="SRR6750041_2_smalltest.fastq"
OUTPUT_DIR="results"
TARGET_MEMORY="8000"
GRANULARITY="100000"
COLLAPSE="true"
ALIGN="star"
KALLISTOINDEXIDX="/mnt/isilon/davidson_lab/ranum/Tools/Kallisto_Index/GRCm38.idx"
KALLISTOINDEXFASTA="/mnt/isilon/davidson_lab/ranum/Tools/Kallisto_Index/Mus_musculus.GRCm38.cdna.all.fa"
STARGENOME="/mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10"
STARGTF="/mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10_Raw/Mus_musculus.GRCm38.93.chr.gtf"

################################
### User Inputs Using Getopt ###
################################
# NOTE on mac systems the options won't work because mac doesnt have the GNU version of getopt by default. GNU getopt can be installed on mac using homebrew. You can do this by running 'brew install gnu_getopt' 
# Once gnu_getopt is installed you can run it with using this '/usr/local/Cellar/gnu-getopt/1.1.6/bin/getopt' as the executable in the place of 'getopt' below.

# read the options
TEMP=`getopt -o n:e:m:1:2:3:f:r:o:t:g:c:a:i:k: --long numcores:,errors:,minreads:,round1barcodes:,round2barcodes:,round3barcodes:,fastqF:,fastqR:,outputdir:,targetMemory:,granularity:,collapseRandomHexamers:,align:,kallistoIndexFasta:,kallistoIndexIDX: -n 'test.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
echo "Checking options..."
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
        -t|--targetMemory)
            case "$2" in
                "") shift 2;;
                *) TARGET_MEMORY=$2 ; shift 2 ;;
            esac ;;
        -g|--granularity)
            case "$2" in
                "") shift 2;;
                *) GRANULARITY=$2 ; shift 2 ;;
            esac ;;
        -c|--collapseRandomHexamers)
            case "$2" in
                "") shift 2;;
                *) COLLAPSE=$2 ; shift 2 ;;
            esac ;;
        -a|--align)
            case "$2" in
                "") shift 2;;
                *) ALIGN=$2 ; shift 2 ;;
            esac ;;
        -i|--kallistoIndexFasta)
            case "$2" in
                "") shift 2;;
                *) KALLISTOINDEXFASTA=$2 ; shift 2 ;;
            esac ;;
        -k|--kallistoIndexIDX)
            case "$2" in
                "") shift 2;;
                *) KALLISTOINDEXIDX=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

###############################
### Write Input Args To Log ###
###############################

if [ $COLLAPSE = true ]
then
ROUND1="Round1_barcodes_new5.txt"
fi

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
echo "targetMemory = $TARGET_MEMORY"
echo "granularity = $GRANULARITY"
echo "collapseRandomHexamers = $COLLAPSE"
echo "align = $ALIGN"
echo "kallistoIndexFasta = $KALLISTOINDEXFASTA"
echo "kallistoIndexIDX = $KALLISTOINDEXIDX"

#######################################
# STEP 1: Demultiplex Using Barcodes  #
#######################################

# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP1: Demultiplex using barcodes. Current time : $now" 

# Demultiplex the fastqr file using barcodes
python demultiplex_using_barcodes.py --minreads $MINREADS --round1barcodes $ROUND1 --round2barcodes $ROUND2 --round3barcodes $ROUND3 --fastqr $FASTQ_R --errors $ERRORS --outputdir $OUTPUT_DIR --targetMemory $TARGET_MEMORY --granularity $GRANULARITY

##########################################################################
# STEP 2: Collapse OligoDT and RandomHexamer Barcodes from the same well #
##########################################################################
# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP2: Collapse OligoDT and RandomHexamer Barcodes from the same well. Current time : $now" 
if [ $COLLAPSE = true ]
then
bash Collapse_RanHex_Odt.sh
fi

##########################################################
# STEP 3: For every cell find matching paired end reads  #
##########################################################
# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP3: Finding read mate pairs. Current time : $now" 

# Now we need to collect the other read pair. To do this we can collect read IDs from the $OUTPUT_DIR files we generated in step one.
# Generate an array of cell filenames
python matepair_finding.py --input $OUTPUT_DIR --fastqf $FASTQ_F --output $OUTPUT_DIR --targetMemory $TARGET_MEMORY --granularity $GRANULARITY

########################
# STEP 4: Extract UMIs #
########################
# Generate a progress message
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "Beginning STEP3: Extracting UMIs. Current time : $now" 

# Implement new method for umi and cell barcode extraction
pushd $OUTPUT_DIR
parallel python3 Extract_BC_UMI.py -F {} -R {}-MATEPAIR ::: $(ls *.fastq)
cat *_1.fastq > MergedCells
parallel rm {} ::: $(*fastq*)
mv MergedCells MergedCells_1.fastq
popd

#rm -r $OUTPUT_DIR-UMI
#mkdir $OUTPUT_DIR-UMI

# Parallelize UMI extraction
#{
#ls $OUTPUT_DIR | grep \.fastq$ | parallel -j $NUMCORES -k "umi_tools extract -I $OUTPUT_DIR/{} --read2-in=$OUTPUT_DIR/{}-MATEPAIR --bc-pattern=NNNNNNNNNN --log=processed.log --read2-out=$OUTPUT_DIR-UMI/{}"
#} &> /dev/null


#################################
# STEP 5: Collect Summary Stats #
#################################
# Print the number of lines and barcode ID for each cell to a file
#echo "$(wc -l results-UMI/*.fastq)" | sed '$d' | sed 's/results-UMI\///g' > linespercell.txt
#Rscript generate_reads_violin.r


###########################
# STEP 6: Perform Mapping #
###########################
# generate batch file
if [ $ALIGN = kallisto ]
then 
    # generate batch file
    rm batch.txt
    for file in $(ls results-UMI/); do
        echo "$(echo $file | sed 's|.fastq||g')" "$(echo $file | sed 's|fastq|umi|g')" "$(echo $file)" >> batch.txt
        python align_kallisto.py -F results-UMI/$file 
    done
    
    pushd results-UMI
    mkdir ../kallisto_output
    kallisto pseudo -i $KALLISTOINDEXIDX -o ../kallisto_output --single --umi -b ../batch.txt
    cp matrix.cells results4_colNames_cellIDs.txt 
    popd

    pushd kallisto_output
    python3 ../prep_TCC_matrix.py -T matrix.tsv -E matrix.ec -O results -I $KALLISTOINDEXFASTA -G geneIDs
    popd
fi

if [ $ALIGN = star ]
then
    pushd $OUTPUT_DIR
    # Run alignment of merged .fastq file using STAR 
    STAR --runThreadN $NUMCORES \
         --readFilesIn MergedCells_1.fastq \
         --outFilterMismatchNoverLmax 0.05 \
         --genomeDir $STARGENOME \
         --alignIntronMax 20000 \
         --outSAMstrandField intronMotif \
         --quantMode GeneCounts \
         --sjdbGTFfile $STARGTF 
 
    # Assign reads to genes
    featureCounts -a $STARGTF \
                  -o gene_assigned \
                  -R BAM Aligned.sortedByCoord.out.bam \
                  -T $NUMCORES

    samtools sort Aligned.sortedByCoord.out.bam -o assigned_sorted.bam
    samtools index assigned_sorted.bam

    # Count UMIs per gene per cell
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv.gz
    popd
fi

#All finished
number_of_cells=$(ls -1 "$OUTPUT_DIR-UMI" | wc -l)
now=$(date '+%Y-%m-%d %H:%M:%S')
echo "a total of $number_of_cells cells were demultiplexed from the input .fastq" 
echo "Current time : $now" 
echo "all finished goodbye" 

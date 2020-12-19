#!/bin/bash

#alias python='python'

#####################################################################################
# Example Use #
# Script modified on 29the august, 2019
# A python script Collapse_Ranhex_Odt.py ( Dipankar / Dumaatravaie ) had replaced the bash script Collapse_Ranhex_Odt.sh
# Fixed the issue with 'Too many arguments error' in parallel ( Dipankar / Dumaatravaie ) in the case of Large Number of cells (> 50,000 or more )
#####################################################################################

#bash splitseqdemultiplex.sh \
# -n 12 \
# -v merged \
# -e 1 \
# -m 10 \
# -1 Round1_barcodes_new4.txt \
# -2 Round2_barcodes_new4.txt \
# -3 Round3_barcodes_new4.txt \
# -f SRR6750041_1_smalltest.fastq \
# -r SRR6750041_2_smalltest.fastq \
# -o results \
# -t 8000 \
# -g 100000 \
# -a star \
# -x /path/to/star/genome/index/folder/GRCm38/ \
# -y /path/to/matching/genome/annotation/gtf/GRCm38.gtf \
# -k /path/to/kallisto/index/.idx/ \
# -i /path/to/kallisto/index/.fasta 

################/media/bachar.d/ec530b5a-02c3-4ebe-8b79-8d8a7fc98220/MirCos_splitSeq/SPLiT-Seq_demultiplexing_annotation_pipeline/New_Script_Test/results_b
# Dependencies #
################
# Python3 must be installed and accessible as "python" from your system's path
type python &>/dev/null || { echo "ERROR python is not installed or is not accessible from the PATH as python"; exit 1; }

# UMI_Tools must be installed and accessible from the PATH as "umi_tools"
type umi_tools &>/dev/null || { echo "ERROR umi_tools is not installed or is not accessible from the PATH as umi_tools"; exit 1; }

# parallel must be installed and accessible from the path as "parallel"
type parallel &>/dev/null || { echo "ERROR parallel is not installed or is not accessible from the PATH as parallel"; exit 1; }

# Not to get the error "/bin/ls: Argument list too long"
# So, a solution is to increase the amount of space available for the stack.
# https://unix.stackexchange.com/questions/45583/argument-list-too-long-how-do-i-deal-with-it-without-changing-my-command
ulimit -s 65536

###########################
### Manually Set Inputs ###
###########################

export NUMCORES="4"
export VERSION="fast"
export ERRORS="1"
export MINREADS="10"
export ROUND1="Round1_barcodes_new5.txt"
export ROUND2="Round2_barcodes_new4.txt"
export ROUND3="Round3_barcodes_new4.txt"
export FASTQ_F="SRR6750041_1_smalltest.fastq"
export FASTQ_R="SRR6750041_2_smalltest.fastq"
export OUTPUT_DIR="results_multiThread"
export TARGET_MEMORY="8000"
export GRANULARITY="100000"
export COLLAPSE="true"
export ALIGN="star"
export STARGENOME="/mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10/"
#STARGTF="GTF /mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10_Raw/Mus_musculus.GRCm38.96.chr.gtf"

export SAF="SAF ../GRCm38_genes.saf"
#KALLISTOINDEXIDX="/mnt/isilon/davidson_lab/ranum/Tools/Kallisto_Index/GRCm38.idx"
#KALLISTOINDEXFASTA="/mnt/isilon/davidson_lab/ranum/Tools/Kallisto_Index/Mus_musculus.GRCm38.cdna.all.fa"

################################
### User Inputs Using Getopt ###
################################
# NOTE on mac systems the options won't work because mac doesnt have the GNU version of getopt by default. GNU getopt can be installed on mac using homebrew. You can do this by running 'brew install gnu_getopt' 
# Once gnu_getopt is installed you can run it with using this '/usr/local/Cellar/gnu-getopt/1.1.6/bin/getopt' as the executable in the place of 'getopt' below.

# read the options
TEMP=`getopt -o n:v:e:m:1:2:3:f:r:o:t:g:c:a:x:y:s:k:i --long numcores:,errors:,minreads:,round1barcodes:,round2barcodes:,round3barcodes:,fastqF:,fastqR:,outputdir:,targetMemory:,granularity:,collapseRandomHexamers:,align:,starGenome:,starGTF:,geneAnnotationSAF:,kallistoIndexIDX:,kallistoIndexFASTA: -n 'test.sh' -- "$@"`
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
        -v|--version)
            case "$2" in
                "") shift 2 ;;
                *) VERSION=$2 ; shift 2 ;;
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
        -x|--starGenome)
            case "$2" in
                "") shift 2;;
                *) STARGENOME=$2 ; shift 2 ;;
            esac ;;
        -y|--starGTF)
            case "$2" in
                "") shift 2;;
                *) STARGTF=$2 ; shift 2 ;;
            esac ;;
        -s|--geneAnnotationSAF)
            case "$2" in
                "") shift 2;;
                *) SAF=$2 ; shift 2 ;;
            esac ;;
        -k|--kallistoIndexIDX)
            case "$2" in
                "") shift 2;;
                *) KALLISTOINDEXIDX=$2 ; shift 2 ;;
            esac ;;
        -i|--kallistoIndexFASTA)
            case "$2" in
                "") shift 2;;
                *) KALLISTOINDEXFASTA=$2 ; shift 2 ;;
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
echo "targetMemory = $TARGET_MEMORY"
echo "granularity = $GRANULARITY"
echo "collapseRandomHexamers = $COLLAPSE"
echo "align = $ALIGN"
echo "starGenome = $STARGENOME"
echo "starGTF = $STARGTF"
echo "geneAnnotationSAF = $SAF"
echo "kallistoIndexIDX = $KALLISTOINDEXIDX"
echo "kallistoIndexFASTA = $KALLISTOINDEXFASTA"

#if [ $COLLAPSE = true ]
#then
#ROUND1="Round1_barcodes_new5.txt"
#fi
if [ $VERSION = fast ]
then
    #######################################
    # STEP 1: Demultiplex Using Barcodes  #
    #######################################
    # Generate a progress message
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP1: Demultiplex using barcodes. Current time : $now" 
    # Demultiplex the fastqr file using barcodes
    mkdir $OUTPUT_DIR
    # Set up a function to parallelize the Demultiplex Using Barcodes Step
    linesInInputFastq=$(wc -l < $FASTQ_R)
    num_linesPerSplitFastq=$(expr $linesInInputFastq / $NUMCORES)
    
    #split --lines=${num_linesPerSplitFastq} $FASTQ_R split_fastq_R_
    #split --lines=${num_linesPerSplitFastq} $FASTQ_F split_fastq_F_
    
    head -n 100 $FASTQ_R > position_learner_fastqr.fastq   
 
    split --number="l/$NUMCORES" $FASTQ_R split_fastq_R_
    split --number="l/$NUMCORES" $FASTQ_F split_fastq_F_

    #my_func() {
    #  python InDevOptimizations/DemultiplexUsingBarcodes_New_V1.py -f "split_fastq_F$1" -r "split_fastq_R$1" -b $GRANULARITY -o $OUTPUT_DIR -e $ERRORS -p -t $MINREADS    
    #  }
    #export -f my_func
    
    #ls split_fastq_F* | awk -F "split_fastq_F" '{print $1}' | parallel my_func {}
    
    ls split_fastq_F* | awk -F "split_fastq_F" '{print $2}' | parallel "python InDevOptimizations/DemultiplexUsingBarcodes_New_V1.py -f split_fastq_F{} -r split_fastq_R{} -b $GRANULARITY -o $OUTPUT_DIR -e $ERRORS -p -t $MINREADS" 
     
    #python InDevOptimizations/DemultiplexUsingBarcodes_New_V1.py -f $FASTQ_F -r $FASTQ_R -b $GRANULARITY -o $OUTPUT_DIR -e $ERRORS -p -t $MINREADS
    #--minreads $MINREADS --round1barcodes $ROUND1 --round2barcodes $ROUND2 --round3barcodes $ROUND3 --fastqr $FASTQ_R --errors $ERRORS --outputdir $OUTPUT_DIR --targetMemory $TARGET_MEMORY --granularity $GRANULARITY
    #echo "$(ls $OUTPUT_DIR/*.fastq | wc -l) results files 'cells'  were demultiplexed from the input .fastq file"
    rm position_learner_fastqr.fastq
    rm split_fastq_*
    ###########################
    # STEP 5: Perform Mapping #
    ###########################
    # generate batch file
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP5: Performing Mapping. Current time : $now" 
    if [ $ALIGN = star ]
    then
        pushd $OUTPUT_DIR
        # Run alignment of merged .fastq file using STAR 

        STAR --runThreadN $NUMCORES \
             --readFilesIn MergedCells_1.fastq \
             --outFilterMismatchNoverLmax 0.05 \
             --genomeDir $STARGENOME \
             --alignIntronMax 20000 \
             --outSAMtype BAM SortedByCoordinate
        
        #cp /media/bachar.d/ec530b5a-02c3-4ebe-8b79-8d8a7fc98220/MirCos_splitSeq/SPLiT-Seq_demultiplexing-master/results/*.bam /media/bachar.d/ec530b5a-02c3-4ebe-8b79-8d8a7fc98220/MirCos_splitSeq/SPLiT-Seq_demultiplexing-master/
    
    if [[ $(echo "$SAF" | awk '{print $1}') = SAF ]]
        then 
        countsMode=$(echo "$SAF" | awk '{print $1}')
        countsFile=$(echo "$SAF" | awk '{print $2}')
    #echo $countsFile
    #echo $countsMode
        # Assign reads to genes
    # Removed -M parameters for excluding multimapping reads
        featureCounts -F $countsMode \
                      -a $countsFile \
                      -o gene_assigned \
                      -R BAM Aligned.sortedByCoord.out.bam \
                      -T $NUMCORES \
                      -M
        else

        countsMode=$(echo "$STARGTF" | awk '{print $1}')
        countsFile=$(echo "$STARGTF" | awk '{print $2}')
    
    #echo "Hello world "
    #echo $countsFile
    #echo "Hello world "
    #echo $countsMode     
    # Removed -M parameters for excluding multimapping reads
    featureCounts -F $countsMode \
                      -a $countsFile \
                      -o gene_assigned \
                      -R BAM Aligned.sortedByCoord.out.bam \
                      -T $NUMCORES \
                      -M
        fi

        samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
        samtools index assigned_sorted.bam

        # Count UMIs per gene per cell
        umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv.gz
        popd
    fi
fi

if [ $VERSION = merged ]
then
    #######################################
    # STEP 1: Demultiplex Using Barcodes  #
    #######################################
    # Generate a progress message
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP1: Demultiplex using barcodes. Current time : $now" 
    # Demultiplex the fastqr file using barcodes
    python demultiplex_using_barcodes.py --minreads $MINREADS --round1barcodes $ROUND1 --round2barcodes $ROUND2 --round3barcodes $ROUND3 --fastqr $FASTQ_R --errors $ERRORS --outputdir $OUTPUT_DIR --targetMemory $TARGET_MEMORY --granularity $GRANULARITY
    echo "$(ls $OUTPUT_DIR/*.fastq | wc -l) results files 'cells'  were demultiplexed from the input .fastq file"

    ##########################################################################
    # STEP 2: Collapse OligoDT and RandomHexamer Barcodes from the same well #
    ##########################################################################
    # Generate a progress message
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP2: Collapse OligoDT and RandomHexamer Barcodes from the same well. Current time : $now" 
    if [ $COLLAPSE = true ]
    then
        # Bash script is painfully slow
	# Replaced by python script, which is 100 times faster then the bash script
	#bash Collapse_RanHex_Odt.sh
	python Collapse_RanHex_Odt.py
    fi
    echo "after collapsing OligoDT and RandomHexamer Barcodes, $(ls $OUTPUT_DIR/*.fastq | wc -l) results files 'cells' remain."

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
    echo "Beginning STEP4: Extracting UMIs. Current time : $now" 

    # Implement new method for umi and cell barcode extraction
    pushd $OUTPUT_DIR
    
    #parallel python ../Extract_BC_UMI.py -R {} -F {}-MATEPAIR ::: $(ls *.fastq)
    # Modified the original parallel command to avoid TOO many arguments error 
    
    ls | grep '\.fastq$' | parallel python ../Extract_BC_UMI.py -R {} -F {}-MATEPAIR
    
    cat *_1.fastq > MergedCells
    
    #parallel rm {} ::: $(ls *fastq*) 
    
    # This command has been changed to avoid Arguments Lists Too Long error
    
    
    # Above command Gives error Argument List Too Long using parallel rm command, and the pipeline halts
    # So, we use either loop to delete the files one by one
    # for i in *.fastq*;do rm "$i";done
    # Or, we can increase stack limit by command ulimit -s 65536, done at the beginning of this script
    # https://unix.stackexchange.com/questions/45583/argument-list-too-long-how-do-i-deal-with-it-without-changing-my-command
    # and we modify the original parallel commands like below, to avoid too many parameters errors
    #parallel rm {} ::: $(ls *fastq*) 
    ls | grep '\.fastq*' | parallel rm {}

    mv MergedCells MergedCells_1.fastq
    popd

    ###########################
    # STEP 5: Perform Mapping #
    ###########################
    # generate batch file
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP5: Performing Mapping. Current time : $now" 
    if [ $ALIGN = star ]
    then
        pushd $OUTPUT_DIR
        # Run alignment of merged .fastq file using STAR 

        STAR --runThreadN $NUMCORES \
             --readFilesIn MergedCells_1.fastq \
             --outFilterMismatchNoverLmax 0.05 \
             --genomeDir $STARGENOME \
             --alignIntronMax 20000 \
             --outSAMtype BAM SortedByCoordinate
        
        #cp /media/bachar.d/ec530b5a-02c3-4ebe-8b79-8d8a7fc98220/MirCos_splitSeq/SPLiT-Seq_demultiplexing-master/results/*.bam /media/bachar.d/ec530b5a-02c3-4ebe-8b79-8d8a7fc98220/MirCos_splitSeq/SPLiT-Seq_demultiplexing-master/
	
	if [ $(echo "$SAF" | awk '{print $1}') = SAF ]
        then 
        countsMode=$(echo "$SAF" | awk '{print $1}')
        countsFile=$(echo "$SAF" | awk '{print $2}')
	#echo $countsFile
	#echo $countsMode
        # Assign reads to genes
	# Removed -M parameters for excluding multimapping reads
        featureCounts -F $countsMode \
                      -a $countsFile \
                      -o gene_assigned \
                      -R BAM Aligned.sortedByCoord.out.bam \
                      -T $NUMCORES \
		      -M
        else

        countsMode=$(echo "$STARGTF" | awk '{print $1}')
        countsFile=$(echo "$STARGTF" | awk '{print $2}')
	
   	#echo "Hello world "
	#echo $countsFile
	#echo "Hello world "
	#echo $countsMode     
	# Removed -M parameters for excluding multimapping reads
	featureCounts -F $countsMode \
                      -a $countsFile \
                      -o gene_assigned \
                      -R BAM Aligned.sortedByCoord.out.bam \
                      -T $NUMCORES \
		      -M
        fi

        samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
        samtools index assigned_sorted.bam

        # Count UMIs per gene per cell
        umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv.gz
        popd
    fi
fi


if [ $VERSION = split ]
then
    #######################################
    # STEP 1: Demultiplex Using Barcodes  #
    #######################################

    # Generate a progress message
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP1: Demultiplex using barcodes. Current time : $now" 

    # Demultiplex the fastqr file using barcodes
    python demultiplex_using_barcodes.py --minreads $MINREADS --round1barcodes $ROUND1 --round2barcodes $ROUND2 --round3barcodes $ROUND3 --fastqr $FASTQ_R --errors $ERRORS --outputdir $OUTPUT_DIR --targetMemory $TARGET_MEMORY --granularity $GRANULARITY


    ##########################################################
    # STEP 2: For every cell find matching paired end reads  #
    ##########################################################
    # Generate a progress message
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "Beginning STEP2: Finding read mate pairs. Current time : $now" 

    # Now we need to collect the other read pair. To do this we can collect read IDs from the $OUTPUT_DIR files we generated in step one.
    # Generate an array of cell filenames
    python matepair_finding.py --input $OUTPUT_DIR --fastqf $FASTQ_F --output $OUTPUT_DIR --targetMemory $TARGET_MEMORY --granularity $GRANULARITY


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
    ls $OUTPUT_DIR | grep \.fastq$ | parallel -j $NUMCORES -k "umi_tools extract -I $OUTPUT_DIR/{} --read2-in=$OUTPUT_DIR/{}-MATEPAIR --bc-pattern=NNNNNNNNNN --log=processed.log --read2-out=$OUTPUT_DIR-UMI/{}"
    } &> /dev/null

    #################################
    # STEP 4: Collect Summary Stats #
    #################################
    # Print the number of lines and barcode ID for each cell to a file
    echo "$(wc -l $OUTPUT_DIR-UMI/*.fastq)" | sed '$d' | sed 's/$OUTPUT_DIR-UMI\///g' > linespercell.txt
    

    ###########################
    # STEP 5: Perform Mapping #
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
        python ../prep_TCC_matrix.py -T matrix.tsv -E matrix.ec -O results -I $KALLISTOINDEXFASTA -G geneIDs
        popd
    fi

number_of_cells=$(ls -1 "$OUTPUT_DIR-UMI" | wc -l)
echo "a total of $number_of_cells cells were demultiplexed from the input .fastq" 
fi


#All finished
#number_of_cells=$(ls -1 "$OUTPUT_DIR-UMI" | wc -l)
now=$(date '+%Y-%m-%d %H:%M:%S')
#echo "a total of $number_of_cells cells were demultiplexed from the input .fastq" 
# Re Initialize the stack limit to default
ulimit -s 8192
echo "Current time : $now" 
echo "all finished goodbye" 


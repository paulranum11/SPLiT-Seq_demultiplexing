#! /bin/bash

# Provide the filenames of the .csv files that contain the barcode sequences. These files should be located in the working directory.
ROUND1="Round1_barcodes_new2.txt"
ROUND2="Round2_barcodes_new2.txt"
ROUND3="Round3_barcodes_new2.txt"


# Provide the filenames of the .fastq files of interest. For this experiment paired end reads are required.
FASTQ_F="SRR6750041_1.fastq"
FASTQ_R="SRR6750041_2.fastq"
FASTQ_TEST="SRR6750041_2_TEST.fastq"


# Add the barcode sequences to a bash array.
declare -a ROUND1_BARCODES=( $(cut -b 1- $ROUND1) )
#printf "%s\n" "${ROUND1_BARCODES[@]}"

declare -a ROUND2_BARCODES=( $(cut -b 1- $ROUND2) )
#printf "%s\n" "${ROUND2_BARCODES[@]}"

declare -a ROUND3_BARCODES=( $(cut -b 1- $ROUND3) )
#printf "%s\n" "${ROUND3_BARCODES[@]}"

# Initialize the counter
count=1

# Make folder for results files
mkdir results


# Search for the barcode in the sample reads file
# Use a for loop to iterate a search for each barcode.  If a match for the first barcode is found search for a match for a second barcode. If a match for the second barcode is found search through the third list of barcodes.
rm ROUND*
rm results/result*



for barcode1 in "${ROUND1_BARCODES[@]}";
    do
    grep -B 1 -A 2 "$barcode1" $FASTQ_R > ROUND1_MATCH.fastq
    echo barcode1.is.$barcode1
    
        if [ -s ROUND1_MATCH.fastq ]
        then

            for barcode2 in "${ROUND2_BARCODES[@]}";
            do
            grep -B 1 -A 2 "$barcode2" ROUND1_MATCH.fastq > ROUND2_MATCH.fastq
            echo barcode2.is.$barcode2
               
                if [ -s ROUND2_MATCH.fastq ]
                then

                    for barcode3 in "${ROUND3_BARCODES[@]}";
                    do
                    grep -B 1 -A 2 "$barcode3" ./ROUND2_MATCH.fastq > ROUND3_MATCH.fastq
                    echo barcode3.is.$barcode3
                    if [ -s ROUND3_MATCH.fastq ]
                    then
                    mv ROUND3_MATCH.fastq results/result.$count.fastq
                    fi

                    count=`expr $count + 1`
                    done
                fi
            done
        fi
    done

find results/ -size  0 -print0 |xargs -0 rm --

#All finished
echo all finished goodbye

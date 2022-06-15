#!/bin/bash
#SBATCH --job-name=AdapterTrim
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    -a GGGGGGGGGG \
#    -g GGAAGCAGTGGTATCAACGCAGAGTGAATGGGAAGCAGTGGTATCAACGC \
#    -o HD_Plate2_5Knuclei_S4_R1_001_adapterTrim.fastq \
#    -p HD_Plate2_5Knuclei_S4_R2_001_adapterTrim.fastq \
#    HD_Plate2_5Knuclei_Fragmented_S4_R1_001.fastq HD_Plate2_5Knuclei_Fragmented_S4_R2_001.fastq

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    --trim-n \
#    -a "A{10}" \
#    -e 0.1 \
#    -o head_fastq_test_trim_R1.fastq \
#    -p head_fastq_test_trim_R2.fastq \
#    head_fastq_test_trim.fastq head_fastq_test_trim.fastq 

cutadapt -f fastq \
	-g N{165} \
	-o output_fromOldVersion/MergedCells_1_trimmed.fastq \
	output_fromOldVersion/MergedCells_1.fastq

#cutadapt -f fastq --match-read-wildcards -m 22 \
#	-g AGGCCAGAGCATTCG \
#	-e 0.1 \
#        -o head_fastq_test_trim_R1.fastq \
#        merged_nanopore_SUP_fastqs_trimmed.fastq 

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    --trim-n \
#    -g "A{10}" \
#    -e 0.1 \
#    -o head_fastq_test_trimIntermediate_R1.fastq \
#     head_fastq_test_trim_R1.fastq 

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    --trim-n \
#    -g "A{10}" \
#    -e 0.1 \
#    -o head_fastq_test_trimIntermediate2_R1.fastq \
#    head_fastq_test_trimIntermediate_R1.fastq

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    --trim-n \
#    -g "T{10}" \
#    -e 0.1 \
#    -o head_fastq_test_trimIntermediate3_R1.fastq \
#    head_fastq_test_trimIntermediate2_R1.fastq

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    --trim-n \
#    -g "T{10}" \
#    -e 0.1 \
#    -o head_fastq_test_trimIntermediate4_R1.fastq \
#    head_fastq_test_trimIntermediate3_R1.fastq

#cutadapt -f fastq --match-read-wildcards -m 22 \
#    --trim-n \
#    -g "T{10}" \
#    -e 0.1 \
#    -o head_fastq_test_trimFinal_R1.fastq \
#    head_fastq_test_trimIntermediate2_R1.fastq

#rm head_fastq_test_trim_R1.fastq 
#rm head_fastq_test_trimIntermediate_R1.fastq
#rm head_fastq_test_trimIntermediate*_R1.fastq

#rm HD_Plate1_5Knuclei_S1_R2_001_adapterTrim.fastq 

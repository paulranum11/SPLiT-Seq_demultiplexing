#!/bin/bash
#SBATCH --job-name=SPLiT_Seq_Demultiplexing
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

python splitseqdemultiplex_0.2.2.py \
	-n 10 \
	-e 1 \
	-m 1 \
	-1 Round1_barcodes_new5.txt \
	-2 Round2_barcodes_new4.txt \
	-3 Round3_barcodes_new4.txt \
	-f SRR6750041_1_smalltest.fastq \
	-r SRR6750041_2_smalltest.fastq \
	-o output \
	-a True \
	-x /mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10 \
	-s GRCm38_genes.saf \
        -b 1000 \
	-p True \
	-l 20000 


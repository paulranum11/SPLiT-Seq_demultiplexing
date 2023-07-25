#!/bin/bash
#SBATCH --job-name=SPLiT_Seq_Demultiplexing
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

python splitseqdemultiplex_0.2.3.py \
	-n 10 \
	-e 1 \
	-m 2 \
	-1 Round1_barcodes_new5.txt \
	-2 Round2_barcodes_new4.txt \
	-3 Round3_barcodes_new4.txt \
	-f /mnt/isilon/davidson_lab/ranum/SequencingAnalysis/HD_S2_AllSamples/HD_S2_1A_S9_R1_001_adapterTrimFinal_head.fastq \
	-r /mnt/isilon/davidson_lab/ranum/SequencingAnalysis/HD_S2_AllSamples/HD_S2_1A_S9_R2_001_adapterTrimFinal_head.fastq \
	-o output \
	-a True \
	-x /mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10 \
	-s GRCm38_genes.saf \
        -i SAF
        -b 10000 \
	-p True \
	-k \
	-l 100000


#1641582384


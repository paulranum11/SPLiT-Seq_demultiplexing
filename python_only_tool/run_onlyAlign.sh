#!/bin/bash
#SBATCH --job-name=SPLiT_Seq_Demultiplexing
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G


python splitseqdemultiplex_onlyAlign.py \
	-t 10 \
	-e 2 \
	-m 10 \
	-1 Round1_barcodes_new5.txt \
	-2 Round2_barcodes_new4.txt \
	-3 Round3_barcodes_new4.txt \
	-f /mnt/isilon/davidson_lab/ranum/SequencingAnalysis/Nanopore_SPLiT-Seq_03_18_2022/HD500Nuclei_merged.fastq \
	-r /mnt/isilon/davidson_lab/ranum/SequencingAnalysis/Nanopore_SPLiT-Seq_03_18_2022/HD500Nuclei_merged.fastq \
	-o output_fromOldVersion \
	-a True \
	-x /mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10 \
	-s GRCm38_genes.saf \
        -b 100000 \
	-p True \
	-l 58414688 \
	-n True \
        -g /mnt/isilon/davidson_lab/ranum/Tools/STAR_Genomes/mm10_Raw/GRCm38.fasta




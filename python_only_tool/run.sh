python splitseqdemultiplex_0.2.2.py \
	-n 15 \
	-e 1 \
	-m 10 \
	-1 Round1_barcodes_new5.txt \
	-2 Round2_barcodes_new4.txt \
	-3 Round3_barcodes_new4.txt \
	-f ~/Desktop/WK14_S2_R1_001_head.fastq \
	-r ~/Desktop/WK14_S2_R2_001_head.fastq \
	-o output \
	-a False \
	-x missing \
	-b 100000 \
	-p True \
	-l 3000000	


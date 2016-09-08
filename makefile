all:
	scripts/quality.pl -f ../transcriptome/CF094_S1_L001_R1_001.fastq -r ../transcriptome/CF094_S1_L001_R2_001.fastq -t 10

clear:
	rm temp/*

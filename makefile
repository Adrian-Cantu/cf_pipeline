all:
	cat makefile
P01:
	cat IDS.txt | xargs -t -I{} ./scripts/quality.pl -f P00_rawreads/{}_R1_001.fastq -r P00_rawreads/{}_R2_001.fastq -t 20 -o {}

clean:
	rm temp/*

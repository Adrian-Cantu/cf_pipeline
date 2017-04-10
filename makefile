all:
	cat makefile
P01:
	cat IDS.txt | xargs -t -I{} ./scripts/quality.pl -f P00_rawreads/{}_R1_001.fastq -r P00_rawreads/{}_R2_001.fastq -t 20 -o {}

P02:
	cat IDS.txt | xargs -I{} sh -c 'pear -f P01_prinseq_output/{}_1.fastq -r P01_prinseq_output/{}_2.fastq -o P02_pear_output/{} -j 20'

P03:
	cat IDS.txt | xargs -t -I{} sh -c 'cat P02_pear_output/{}.assembled.fastq P01_prinseq_output/{}_1_singletons.fastq P01_prinseq_output/{}_2_singletons.fastq P02_pear_output/{}.unassembled.forward.fastq | paste - - - - | cut -f 1,2 | sed "s/^@/>/" | tr "\t" "\n" > P03_qfpassed_fasta/qc_{}.fasta'

P04:
	cat IDS.txt  | xargs -I{} sh -c '/home/vito/TOOLS/bin/smalt-0.7.6/src/smalt map -y 0.9  -n 20 -f samsoft -o P04_map_HG/mapHG_{}.samsoft /home/vito/TOOLS/DBS/hg19 P03_qfpassed_fasta/qc_{}.fasta'

P045:
	cat IDS.txt | xargs -t -I{} sh -c 'perl /home/vito/TOOLS/bin/samsoft2fasta_nohits.pl P04_map_HG/mapHG_{}.samsoft > P04_noHG/noHG_qf_{}.fasta'


P05:
	cat IDS.txt | xargs -I{} sh -c '/home/vito/TOOLS/bin/smalt-0.7.6/src/smalt map -y 0.8 -n 20 -f samsoft -o P05_map_rRNA/map_SILVA_{}.samsoft /home/vito/TOOLS/DBS/SILVA_all P04_noHG/noHG_qf_{}.fasta'

P055:
	cat IDS.txt | xargs -t -I{} sh -c 'perl /home/vito/TOOLS/bin/samsoft2fasta_nohits.pl P05_map_rRNA/map_SILVA_{}.samsoft > P05_noHG_norRNA/norRNA_noHG_qf_{}.fasta'

P07:
	cat IDS.txt | xargs -I{} sh -c '/home/vito/TOOLS/bin/smalt-0.7.6/src/smalt map -y 0.8 -n 20 -f samsoft -o P07_map_viral_revseq/map_viralrefseq_{}.samsoft /home/vito/TOOLS/DBS/viralrefseq  P05_noHG_norRNA/norRNA_noHG_qf_{}.fasta'
clean:
	rm temp/*

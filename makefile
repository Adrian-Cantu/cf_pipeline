all:
	cat makefile
P01:
	cat IDS3.txt IDS4_bis.txt IDS6.txt | xargs -t -I{} ./scripts/quality.pl -f P00_rawreads/{}_R1_001.fastq -r P00_rawreads/{}_R2_001.fastq -t 8 -o {}


P015: 
	cat ${IDS} | xargs -i -t sh -c  'hisat2 -x ~/db/hisat_db/grch38/genome -1  P01_prinseq_output/{}_1.fasta -2  P01_prinseq_output/{}_2.fasta -S P015_hisat_output/{}_hg.sam -f  --new-summary --time --summary-file  P015_hisat_output/{}_log'



P016: P015
	cat ${IDS} |  xargs -i -t sh -c 'samtools view -S -f  76  P015_hisat_output/{}_hg.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\t" "\n" > P016_hisat_nohit/{}_no_hit_R1.fasta ; samtools view -S -f  140  P015_hisat_output/{}_hg.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\t" "\n" > P016_hisat_nohit/{}_no_hit_R2.fasta'

P017: P016
	cat ${IDS} | xargs -i -t sh -c  'hisat2 -x ~/db/univec/univec -1  P016_hisat_nohit/{}_no_hit_R1.fasta -2  P016_hisat_nohit/{}_no_hit_R2.fasta -S P017_hisat_univec_output/{}_univec.sam -f  --new-summary --time --summary-file  P017_hisat_univec_output/{}_log'


P018: P017
	cat ${IDS} |  xargs -i -t sh -c 'samtools view -S -f  76  P017_hisat_univec_output/{}_univec.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\t" "\n" > P018_hisat_univec_nohit/{}_polish_R1.fasta ; samtools view -S -f  140  P017_hisat_univec_output/{}_univec.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\t" "\n" > P018_hisat_univec_nohit/{}_polish_R2.fasta'


P019: P018
	cat ${IDS} | xargs -i -t sh -c  'hisat2 -x ~/db/viral_refseq/viral_refseq -U  P018_hisat_univec_nohit/{}_polish_R1.fasta -S P019_hisat_viral_refseq/{}_viral_refseq.sam -f  --new-summary --time --summary-file  P019_hisat_viral_refseq/{}_log'


P020: P019
	cat ${IDS} | xargs -t -i sh -c 'grep -v ^@ P019_hisat_viral_refseq/{}_viral_refseq.sam | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\t"  > P020_viral_hits/{}_hits_viral_refseq.tab '

table: P020
	cat ${IDS} |  xargs -i sh -c 'echo -n "{} " ; grep -c ">" P018_hisat_univec_nohit/{}_polish_R1.fasta ' > Tj.txt
P021: table 
	perl frap_normalization.pl -t Tj.txt -m -f ~/db/viral_refseq/viral_refseq.fasta  P020_viral_hits/*.tab > P021_table_viral_hits/frap_viral_refseq.tab

P023:
	cat ${IDS} |  xargs -i sh -c 'python ~/SUPER-FOCUS/superfocus.py -q P018_hisat_univec_nohit/{}_polish_R1.fasta -dir P023_super_focus -a diamond -o {}'

focus:
	python ~/SUPER-FOCUS/superfocus.py -m 1 -q focus_link -dir P023_super_focus -a diamond -o all



test:
	echo ${IDS}


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
P075:
	cat IDS.txt | xargs -t -I{} sh -c 'cat P07_map_viral_refseq/map_viralrefseq_{}.samsoft | grep -v ^@ | cut -f 3 | sort -n | uniq -c | sort -nr |  sed -e "'"s/^ *//"'" | tr "'" "'" \\t | paste -- >  P07_map_viral_refseq/hits_viralrefseq_{}.tab'

P12:
	cat IDS.txt | xargs -I{} sh -c '/home1/acantu/share/smalt map -y 0.8 -n 20 -f samsoft -o P13_spades_map/map_all_cat_{}.samsoft P12_spades_cross/ALL_CONTIGS  P05_noHG_norRNA/norRNA_noHG_qf_{}.fasta'

P13:
	cat IDS.txt | xargs -t -I{} sh -c 'cat P13_spades_map/map_all_cat_{}.samsoft | grep -v ^@ | cut -f 3 | sort -n | uniq -c | sort -nr |  sed -e "'"s/^ *//"'" | tr "'" "'" \\t | paste -- >  P13_spades_map/hits_contigs_{}.tab'

clean:
	rm temp/*

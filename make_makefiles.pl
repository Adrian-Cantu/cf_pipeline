#!/usr/bin/perl
#
use strict;
use warnings;
my $header = <<"EOF";
hisat2 = /home1/acantu/share/hisat2-2.1.0/hisat2
db_folder = /home1/acantu/share/db
focus_folder = /home1/acantu/bin/SUPERFOCUS_0.27


EOF

open INFILE, '<', $ARGV[0] || die $!;;
my @order;

while(<INFILE>) {
	chomp;
	push @order,$_;
}
close(INFILE);

open OUT, '>', "P01_prinseq_output/makefile" || die $!;
print OUT $header;
#open OUT2, '>', "P015_hisat_output/makefile" || die $!;
#print OUT2 $header;
#open OUT3, '>', "" || die $!;
#print OUT3 $header;

open MAIN , '>' , "kk_makefile" || die $!;
print MAIN $header;
my %hash_single;
my %hash_double;
foreach(@order) {
	my @f=split;
	if (scalar(@f)==2) {
		print OUT "$f[0].fasta : ../P00_rawreads/$f[1]";
		print OUT "\tperl ../prinseq-lite.pl -verbose -fastq ../P00_rawreads/$f[1] -derep 1245 -lc_method entropy -lc_threshold 50 -trim_qual_right 15 -trim_qual_left 15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -trim_tail_left 5 -trim_tail_right 5 -min_len 60 -min_qual_mean 20 -ns_max_p 1 -rm_header  -out_bad null -out_format 1 -out_good ../P01_prinseq_output/$f[0]\n\n";
	} elsif (scalar(@f)==3) {
		print MAIN "P01_prinseq_output/$f[0]_1%fasta P01_prinseq_output/$f[0]_2%fasta P01_prinseq_output/$f[0]_1_singletons%fasta P01_prinseq_output/$f[0]_2_singletons%fasta : P00_rawreads/$f[1] P00_rawreads/$f[2]\n";
		print MAIN "\tscripts/quality\$*pl -f P00_rawreads/$f[1] -r P00_rawreads/$f[2] -t 8 -o $f[0]\n\n";
		print MAIN "P015_hisat_output/$f[0]_hg.sam : P01_prinseq_output/$f[0]_1.fasta P01_prinseq_output/$f[0]_2.fasta\n";
		print MAIN "\t\${hisat2} -x \${db_folder}/grch38/genome -1 P01_prinseq_output/$f[0]_1.fasta -2 P01_prinseq_output/$f[0]_2.fasta -S P015_hisat_output/$f[0]_hg.sam -f  --new-summary --time --summary-file  P015_hisat_output/$f[0]_log\n\n";

        print MAIN "P016_hisat_nohit/$f[0]_no_hit_R1%fasta P016_hisat_nohit/$f[0]_no_hit_R2%fasta : P015_hisat_output/$f[0]_hg.sam\n";
        print MAIN qq(\tsamtools view -S -f  76  P015_hisat_output/$f[0]_hg.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P016_hisat_nohit/$f[0]_no_hit_R1.fasta\n);
        print MAIN qq(\tsamtools view -S -f  140  P015_hisat_output/$f[0]_hg.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P016_hisat_nohit/$f[0]_no_hit_R2.fasta\n\n);

        print MAIN "P017_hisat_univec_output/$f[0]_univec.sam  : P016_hisat_nohit/$f[0]_no_hit_R1.fasta P016_hisat_nohit/$f[0]_no_hit_R2.fasta\n";
        print MAIN "\t\${hisat2} -x \${db_folder}/univec/univec -1  P016_hisat_nohit/$f[0]_no_hit_R1.fasta -2  P016_hisat_nohit/$f[0]_no_hit_R2.fasta -S P017_hisat_univec_output/$f[0]_univec.sam -f  --new-summary --time --summary-file  P017_hisat_univec_output/{}_log\n\n";

        print MAIN "P018_hisat_univec_nohit/$f[0]_polish_R1%fasta  P018_hisat_univec_nohit/$f[0]_polish_R2%fasta : P017_hisat_univec_output/$f[0]_univec.sam\n";
        print MAIN qq(\tsamtools view -S -f  76  P017_hisat_univec_output/$f[0]_univec.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P018_hisat_univec_nohit/$f[0]_polish_R1.fasta\n);
        print MAIN qq(\tsamtools view -S -f  140  P017_hisat_univec_output/$f[0]_univec.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P018_hisat_univec_nohit/$f[0]_polish_R2.fasta\n\n);

        print MAIN "P019_hisat_viral_refseq/$f[0]_viral_refseq.sam : P018_hisat_univec_nohit/$f[0]_polish_R1.fasta  P018_hisat_univec_nohit/$f[0]_polish_R2.fasta\n"; 
        print MAIN "\t\${hisat2} -x \${db_folder}/viral_refseq -U  P018_hisat_univec_nohit/$f[0]_polish_R1.fasta -S P019_hisat_viral_refseq/$f[0]_viral_refseq.sam -f --new-summary --time --summary-file P019_hisat_viral_refseq/$f[0]_log\n\n";

        print MAIN "P020_viral_hits/$f[0]_hits_viral_refseq.tab : P019_hisat_viral_refseq/$f[0]_viral_refseq.sam\n";
        print MAIN qq(\tgrep -v ^@ P019_hisat_viral_refseq/$f[0]_viral_refseq.sam | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\\t"  > P020_viral_hits/$f[0]_hits_viral_refseq.tab\n\n);
	}
}


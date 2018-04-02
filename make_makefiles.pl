#!/usr/bin/perl
#
use strict;
use warnings;
my $header = <<"EOF";
hisat2 = /home1/acantu/share/hisat2-2.1.0/hisat2
db_folder = /home1/acantu/share/db
focus_folder = /home1/acantu/bin/SUPERFOCUS_0.27
threads = 20

all : result_$ARGV[0]
	\@echo "DONE :)" 

makefile: $ARGV[0]
	./make_makefiles.pl $ARGV[0]

EOF

open INFILE, '<', $ARGV[0] || die $!;;
my @order;

while(<INFILE>) {
	chomp;
	push @order,$_;
}
close(INFILE);

open MAIN , '>' , "makefile" || die $!;
print MAIN $header;
select(MAIN);
my @single_id;
my @double_id;
foreach(@order) {
	my @f=split;
	if (scalar(@f)==2) {
		print "P01_prinseq_output/$f[0]_good_out.fasta : P00_rawreads/$f[1]\n";
		print "\tprinseq++ -fastq P00_rawreads/$f[1] -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 15  -rm_header  -out_name P01_prinseq_output/$f[0] -threads 40 -out_format 1\n\n";
	    print "P02_map_human/$f[0]_hg.sam: P01_prinseq_output/$f[0]_good_out.fasta\n";
        print "\t\${hisat2} -x \${db_folder}/grch38/genome -p \${threads} -U P01_prinseq_output/$f[0]_good_out.fasta -S P02_map_human/$f[0]_hg.sam -f  --new-summary --time --summary-file  P02_map_human/$f[0]_log\n\n";
        print "P02_map_human/$f[0]_no_hit.fasta : P02_map_human/$f[0]_hg.sam\n";
        print qq(\tsamtools view -S -f  4  P02_map_human/$f[0]_hg.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P02_map_human/$f[0]_no_hit.fasta\n);

        print "P03_map_univec/$f[0]_univec.sam : P02_map_human/$f[0]_no_hit.fasta\n";
        print "\t\${hisat2} -x \${db_folder}/univec/univec -p \${threads} -U  P02_map_human/$f[0]_no_hit.fasta -S P03_map_univec/$f[0]_univec.sam -f  --new-summary --time --summary-file  P03_map_univec/$f[0]_log\n\n";

        print "P03_map_univec/$f[0]_no_hit.fasta: P03_map_univec/$f[0]_univec.sam\n";
        print qq(\tsamtools view -S -f  4  P03_map_univec/$f[0]_univec.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P03_map_univec/$f[0]_no_hit.fasta\n\n);

        print "P04_map_gg/$f[0]_gg.sam :  P03_map_univec/$f[0]_no_hit.fasta\n";
        print "\t\${hisat2} -x \${db_folder}/greengenes/gg -p \${threads}  -U  P03_map_univec/$f[0]_no_hit.fasta -S P04_map_gg/$f[0]_gg.sam  -f  --new-summary --time --summary-file  P03_map_univec/$f[0]_log\n\n";

        print "P05_polish/$f[0]_polish.fasta : P04_map_gg/$f[0]_gg.sam\n";
        print qq(\tsamtools view -S -f  4  P04_map_gg/$f[0]_gg.sam | cut -f 1,10 | sort | sed "s/^/>/" | tr "\\t" "\\n" > P05_polish/$f[0]_polish.fasta\n\n);
  #      print "P05_polish/$f[0]_polish.fasta  : P04_map_gg/$f[0]_no_hit.fasta"
        print MAIN "P11_virus/$f[0]_viral_refseq.sam : P05_polish/$f[0]_polish.fasta\n";
        print MAIN "\t\${hisat2} -x \${db_folder}/viral_refseq -p \${threads} -U  P05_polish/$f[0]_polish.fasta -S P11_virus/$f[0]_viral_refseq.sam -f --new-summary --time --summary-file P11_virus/$f[0]_log\n\n";
        print MAIN "P11_virus/$f[0]_hits_viral_refseq.tab : P11_virus/$f[0]_viral_refseq.sam\n";
        print MAIN qq(\tgrep -v ^@ P11_virus/$f[0]_viral_refseq.sam | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\\t"  > P11_virus/$f[0]_hits_viral_refseq.tab\n\n);
        push @double_id,$f[0];
        

##################################################################################################   
    } elsif (scalar(@f)==3) {
		print MAIN "P01_prinseq_output/$f[0]_good_out_R1%fasta P01_prinseq_output/$f[0]_good_out_R2%fasta : P00_rawreads/$f[1] P00_rawreads/$f[2]\n";
#		print MAIN "\tscripts/quality\$*pl -f P00_rawreads/$f[1] -r P00_rawreads/$f[2] -t 8 -o $f[0]\n\n";
        print MAIN "\tprinseq++ -fastq P00_rawreads/$f[1] -fastq2 P00_rawreads/$f[2] -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 15  -rm_header  -out_name P01_prinseq_output/$f[0] -threads 40 -out_format 1  \n\n";
		print MAIN "P015_hisat_output/$f[0]_hg.sam : P01_prinseq_output/$f[0]_good_out_R1.fasta P01_prinseq_output/$f[0]_good_out_R2.fasta\n";
		print MAIN "\t\${hisat2} -x \${db_folder}/grch38/genome -1 P01_prinseq_output/$f[0]_good_out_R1.fasta -2 P01_prinseq_output/$f[0]_good_out_R2.fasta -S P015_hisat_output/$f[0]_hg.sam -f  --new-summary --time --summary-file  P015_hisat_output/$f[0]_log\n\n";

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
        push @double_id,$f[0];

	}
}
    my @tab_list= map {'P11_virus/'.$_.'_hits_viral_refseq.tab'} @double_id ;
    my $tab_list= join(' ',@tab_list);
    print MAIN "Tj.txt : $tab_list\n";
    print MAIN qq(\tcat $ARGV[0] | cut -f 1 | xargs -I{} sh -c 'echo -n "{} " ; grep -c ">" P018_hisat_univec_nohit/{}_polish_R1.fasta ' > Tj.txt\n\n);

    print MAIN "result_$ARGV[0] : Tj.txt\n";
    print MAIN qq(\tperl frap_normalization.pl -t Tj.txt -m -f \${db_folder}/viral_refseq.fasta $tab_list > result_$ARGV[0]);



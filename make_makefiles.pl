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

my %hash_single;
my %hash_double;
foreach(@order) {
	my @f=split;
	if (scalar(@f)==2) {
	} elsif (scalar(@f)==3) {
		(my $file1 = $f[1]) =~ s/\./%/;
		(my $file2 = $f[2]) =~ s/\./%/;
		print OUT "$f[0]_1%fasta $f[0]_2%fasta $f[0]_1_singletons%fasta $f[0]_2_singletons%fasta : ../P00_rawreads/$f[1] ../P00_rawreads/$f[2]\n";
		print OUT "\t../scripts/quality\$2*pl -f ../P00_rawreads/$f[1] -r ../P00_rawreads/$f[2] -t 8 -o $f[0]\n\n";
	}
}


#! /usr/bin/env perl
#
# Short description for P07_make_viral_refseq_table.pl
#
# Author vito <vito@rambox>
# Version 0.1
# Copyright (C) 2017 vito <vito@rambox>
# Modified On 2017-04-10 10:26
# Created  2017-04-10 10:26
#
use strict;
use warnings;

my @glob_ARGV=@ARGV;
my @order;
my %sample_id;
my %hits;
my $gene;


foreach my $file_name (@glob_ARGV) {
	open IN, '<', $file_name;
	#my $sample=<IN>;
	#chomp($sample);
	$sample_id{$file_name}=(split(/\_/,$file_name))[-1];
	$hits{$file_name}={};
	while(<IN>) {
		chomp;
		my @fields=split;
		push @order,$fields[1] unless $fields[1]~~@order;
		$hits{$file_name}->{$fields[1]}=$fields[0];
	}
	close(IN);
}

print "-";
foreach my $file_name (@glob_ARGV) {
	print "\t$sample_id{$file_name}"
}
print "\n";


foreach(@order) {
	next if $_ eq '*';
	$gene = $_;
	print "$gene";
	foreach(@glob_ARGV) {
		if (exists $hits{$_}->{$gene}) {	
			print "\t$hits{$_}->{$gene}";
		} else {
			print "\t0";
		}

	}		
	print "\n";
}





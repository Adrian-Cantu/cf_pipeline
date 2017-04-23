#! /usr/bin/env perl
#
# Short description for translate_viralDB.pl
#
# Author vito <vito@rambox>
# Version 0.1
# Copyright (C) 2017 vito <vito@rambox>
# Modified On 2017-04-10 13:17
# Created  2017-04-10 13:17
#
use strict;
use warnings;
my %vir_trans;

open DB,'<',"translate_viralDB.tsv";
open IN,'<',$ARGV[0];

while (<DB>) {
	chomp;
	my @field=split;
	my $id = shift(@field);
	my $desc = join(' ',@field);
	$vir_trans{$id}=$desc;
}
my $line=<IN>;
print $line;

while (<IN>) {
	my @field=split;
	my $id = shift(@field);
	my $rest = join("\t",@field);
	print "'$vir_trans{$id}'\t$rest\n";
}




#!/usr/bin/perl
#
use strict;
use warnings;

my %hash;

my $list=`ls -1p P00_rawreads | grep fastq\$`;
chomp $list;
my @list=split(/\n/,$list);
foreach (@list) {
	my @f=split(/_/,$_);
	my $id=join('_',$f[0],$f[1],$f[2]);
	if (exists $hash{$id}) {
		push @{$hash{$id}}, $_;
	} else {
		$hash{$id}=[$_];
	}
}

foreach (keys %hash) {
	print join("\t",$_,@{$hash{$_}}) , "\n";
}

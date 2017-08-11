#!/usr/bin/perl
#
use Getopt::Long;
use strict;
use warnings;

my $id;
my $hits_only;
my $norm_only;
my %desc;
my %Tj;
my $tj_file;
my $millon;
GetOptions ("f=s"     => \$id,
	    "h"       => \$hits_only,	
	    "n"       => \$norm_only,    
	    "t=s"       => \$tj_file,
 	    "m"       => \$millon,	    
);  # flag

my @glob_ARGV=@ARGV;
my %list;
my %length;
my %sum;
my %norm;
my %samples;
my @order;
open IDS , '<' , $id;

#start hashs
foreach (@ARGV){
	$list{$_}={};
	$norm{$_}={};
	$sum{$_}=0;
}

#read the fasta and extract id, length and description for each sequence
my $ll;
my $ll_name;
while(<IDS>) {
	if (/^>/) {
		chomp;
		my @f=split;
		my $id = shift @f;
		$id =~ s/^>//;
		$desc{$id}=join(' ',@f);
		push @order,$id;

		$length{$ll_name}=length($ll) if $ll_name;
		$ll='';
		$ll_name=$id;
		foreach (@ARGV){
        		$list{$_}->{$id}=0;
		}
	}
	else {
		chomp;
		$ll=$ll.$_;
	}
}
$length{$ll_name}=length($ll);


my %filehandlers;
foreach(@ARGV) {
	chomp;
	open $filehandlers{$_}, '<' , $_  or die "Can't open $_ for output: $!";
	my $filename = (split(/\//,$_))[-1];
	my $sample_id = (split(/\./,$filename))[0];
	my @sample = split(/\_/,$sample_id);
	
	#getting IDS from filename, this need to be much more general

	my $kk = pop @sample;
	$kk = pop @sample;
	$kk = pop @sample;
	$kk = join('_',@sample);
	$samples{$_} = $kk;	
}


open TJ , '<', $tj_file;
while(<TJ>) {
	chomp;
	my @f=split;
	$Tj{$f[0]}=$f[1];
}


foreach(@ARGV) {
	my $f_id=$_;
	while(readline($filehandlers{$f_id})) {
#next if $_=~/^\-$/;
	my @fields=split;
	print "$fields[1]\n" unless exists $list{$f_id}->{$fields[1]}; 
	$list{$f_id}->{$fields[1]}=$list{$f_id}->{$fields[1]}+$fields[0];
	}
}
foreach(@glob_ARGV) {
	my $f_id=$_;
	my $total=0;
	foreach(@order) {
		$total=$total+$list{$f_id}->{$_};
	}
	$sum{$f_id}=$total;
}
		
foreach(@glob_ARGV) {
	my $f_id=$_;
	foreach(@order) {
		$norm{$f_id}->{$_}=($list{$f_id}->{$_}/$Tj{$samples{$f_id}})*(300/$length{$_});
	}
}

print "-";
foreach(@glob_ARGV) {
	print "\t$samples{$_}"
}
print "\n";

foreach(@order) {
	my $g_id=$_;
	print "$g_id";
	foreach(@glob_ARGV) {
        	my $f_id=$_;
		if ($hits_only) {
			print "\t$list{$_}->{$g_id}";
		} elsif ($norm_only) {
			print "\t$norm{$f_id}->{$g_id}";
		} elsif ($millon) {
			my $t_kk = $norm{$f_id}->{$g_id} * 1000000;
		       print "\t$t_kk";
	       }	       
		
	}
	print "\t$desc{$g_id}\n";
}


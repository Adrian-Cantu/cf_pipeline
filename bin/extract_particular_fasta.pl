#!/usr/bin/perl -w

use strict;
my $p=0;
my $a;
my $line;
my $file = shift @ARGV;
open IN, '<' , $file;

while(<IN>) {
    $p=0 if $_=~/^\>/;
    $line=$_;
    foreach $a (@ARGV) {
        $p=1 if $line=~/$a/;
    }
    print $_ if $p;
} 




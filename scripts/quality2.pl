#! /usr/bin/env perl
#
# Short description for quality.pl
#
# Author Vito Adrian Cantu <garbanyo@gmail.com>
# Version 0.1
# Copyright (C) 2016 Vito Adrian Cantu <garbanyo@gmail.com>
# Modified On 2016-08-24 13:24
# Created  2016-08-24 13:24
#
use strict;
use warnings;
use Getopt::Long;
use 5.010;

my $help;
my $forward_file;
my $reverse_file;
my $lines_R;
my $lines_F;
my $threads;
my $out_file_prefix;
GetOptions (
    'help'         => \$help,
    'reverse:s'    => \$reverse_file, 
    'forward=s'    => \$forward_file,
    'threads=i'    => \$threads,
    'out=s'        => \$out_file_prefix,	
);
&help if ($help);

$threads=int($threads);
$lines_F = `echo \$((\$((\$((\$((\$(wc -l < $forward_file) / 4  )) / $threads )) + 1 )) * 4)) `;
if ($reverse_file ne "") {
    $lines_R = `echo \$((\$((\$((\$((\$(wc -l < $reverse_file) / 4  )) / $threads )) + 1 )) * 4)) `;
}

print "lines f $lines_F \nlines r $lines_R\n";
chomp $lines_F;
chomp $lines_R;


my @chars = ("A".."Z", "a".."z");
my $string;
$string .= $chars[rand @chars] for 1..4;
print "writing to $string";
#system $run;


mkdir "../temp" unless -d "../temp";
my $run= "split -a 3 -l $lines_F -d $forward_file ../temp/$string" . "_R1.fastq";
system $run;

if ($reverse_file ne "") {
    $run= "split -a 3 -l $lines_R -d $reverse_file ../temp/$string" . "_R2.fastq";
    system $run;
    }

#  $run="perl prinseq-lite.pl -verbose -fastq $1 -fastq2 $2 -derep 1245 -lc_method entropy -lc_threshold 50 -trim_qual_right 15 -trim_qual_left 15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_rule lt -trim_qual_window 2 -trim_tail_left 5 -trim_tail_right 5 -min_len 60 -min_qual_mean 15 -ns_max_p 1 -rm_header  -out_bad null -out_good 'passed'

my $n = $threads;
my $forks = 0;
for (1 .. $n) {
  my $pid = fork;
  if (not defined $pid) {
     warn 'Could not fork';
     next;
  }
  if ($pid) {
    $forks++;
    printf "In the parent process PID ($$), Child pid: $pid Num of fork child processes: $forks\n";
  } else {
    printf "In the child process PID ($$)\n"; 
    sleep 2;
    my $f_num = sprintf("%03d", $forks);
    system "perl ../prinseq-lite.pl -verbose -fastq ../temp/$string" . "_R1.fastq$f_num -fastq2 ../temp/$string" . "_R2.fastq$f_num -derep 1245 -lc_method entropy -lc_threshold 50 -trim_qual_right 15 -trim_qual_left 15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -trim_tail_left 5 -trim_tail_right 5 -min_len 60 -min_qual_mean 20 -ns_max_p 1 -rm_header  -out_bad null -out_format 1 -out_good '../temp/passed_$string" . "$f_num' 2> ../temp/log_$string" . "$f_num";
    say "Child ($$) exiting";
    exit;
  }
}
 
for (1 .. $forks) {
   my $pid = wait();
   say "Parent saw $pid exiting";
}

system "cat ../temp/passed_$string*_1.fasta > $out_file_prefix"."_1.fasta";
system "cat ../temp/passed_$string*_2.fasta > $out_file_prefix"."_2.fasta";
system "cat ../temp/passed_$string*_1_singletons.fasta > $out_file_prefix"."_1_singletons.fasta";
system "cat ../temp/passed_$string*_2_singletons.fasta > $out_file_prefix"."_2_singletons.fasta";
#system "make clean";



sub help{
    my $helpT = qq(


-h print this help page

for more information use "perldoc <quality.pl>"
);
    print $helpT;
    exit;
    }

=head1 VERSION

This documentation refers to <quality.pl> version 0.1

=head1 SYNOPSIS

-

=head1 USE

-

=head1 OPTIONS

-

=head1 DESCRIPTION

-

=head1 CONFIGURATION AND ENVIRONMENT

-

=head1 DEPENDENCIES

-

=head1 INCOMPATIBILITIES

-

=head1 BUGS AND LIMITATIONS

There are no known bugs in this app.
Please report problems to Vito Adrian Cantu <garbanyo@gmail.com>

=head1 AUTHOR

Vito Adrian Cantu <garbanyo@gmail.com>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2016 Vito Adrian Cantu <garbanyo@gmail.com>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
~


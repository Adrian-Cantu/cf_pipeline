#!/bin/bash

perl prinseq-lite.pl -verbose -fastq $1 -fastq2 $2 -derep 1245 -lc_method entropy -lc_threshold 50 -trim_qual_right 15 -trim_qual_left 15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_rule lt -trim_qual_window 2 -trim_tail_left 5 -trim_tail_right 5 -min_len 60 -min_qual_mean 15 -ns_max_p 1 -rm_header  -out_bad null -out_good 'passed'

pear -f passed_1.fastq -r passed_2.fastq  -o  pair_reads  -j 20


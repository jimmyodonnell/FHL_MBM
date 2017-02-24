#!/usr/bin/env bash

the_file="/Users/threeprime/Dropbox/FHL_MBM/Blast_iterative/Blast_Preclusters_GenBank/FHL_blast_GenBank_results.txt"

outfile="/Users/threeprime/Documents/Projects/FHL_MBM/data/blast/genbank/blast_results_nogit.txt"

awk -F'\t' 'NF == 14 { print $0 }' $the_file > temp_nf.tmp

awk -F'\t' ' { print $3 }' temp_nf.tmp | \
  LC_ALL='C' sort | \
  LC_ALL='C' uniq | \
  LC_ALL='C' grep ^[^0-9] \
  > temp_badchar.tmp

awk -F'\t' ' { print $11 }' temp_nf.tmp | \
  LC_ALL='C' sort | \
  LC_ALL='C' uniq | \
  LC_ALL='C' grep ^[^0-9] \
  >> temp_badchar.tmp

LC_ALL='C' fgrep -v -f temp_badchar.tmp temp_nf.tmp > $outfile

rm temp_nf.tmp


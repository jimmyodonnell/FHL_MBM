#!/usr/bin/env bash

# strip taxon -- this could always change
awk -F'|' '{ 
  if($0 ~ /^>/) 
    print $1; 
  else 
    print $0 
}' FHL_plates1-3_aligned.fasta > FHL_plates_1-3_aligned_notax.fasta

# combine files
cat FHL_plates_1-3_aligned_notax.fasta FHL_Other_plates.fasta > combined.fasta

# strip indels ('-')
awk '{
  if(!/^>/) 
    gsub("-",""); 
  print $0
}' combined.fasta > FHL_mbm_co1.fasta

# rm FHL_plates_1-3_aligned_notax.fasta
# rm combined.fasta

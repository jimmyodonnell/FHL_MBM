#!/usr/bin/env bash

# 
SEQ_TAX_CSV="data/fhl_co1/fhl_mbm_co1_taxa.csv"

NEW_FASTA="data/fhl_co1/FHL_Other_plates.fasta"

awk 'BEGIN {OFS="";} /^>/ { 
  gsub(">","") ; print $0,"," 
}' "${NEW_FASTA}" >> "${SEQ_TAX_CSV}"
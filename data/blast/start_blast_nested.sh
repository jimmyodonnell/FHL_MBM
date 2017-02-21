#!/usr/bin/env sh
module load tools/local

my_dir="/pool/genomics/leraym/FHL_Blast_iterative/Fastas"
for infile in ${my_dir}/*.fasta; do
   echo $infile
   # Uses the script q-wait to only submit more jobs when <njobs jobs are running
   # -njobs: max number of jobs to run at once
   # "FHL_blast_test": text job name to filter by
   q-wait -njobs 500 FHL_blast_test 
   qsub blast_nested.job "${infile}" "1e-50 1e-45 1e-40 1e-35 1e-30 1e-25 1e-20"
done

date
echo "All jobs have now been submitted"

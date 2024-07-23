#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/ASEReadCounter_missense-new_batch/

for SAMPLE in ${INPUT_DIR}/*_ASE.output.table
  do 

  echo " BEGIN COUNTS TO SNPS FOR: $SAMPLE"
  bsub -J countsToSnps -P acc_ISDS -q premium -n 4 -R rusage[mem=16000] -M 80000 -W 10:00 -o logs/countsToSnps_%J.stdout -e logs/countsToSnps_%J.stderr ./4_run_counts_to_snps.sh ${SAMPLE}

done

## EOF ##
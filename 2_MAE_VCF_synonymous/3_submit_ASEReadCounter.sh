#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_DUP_REMOVED_synonymous-new_batch

for SAMPLE in ${INPUT_DIR}/*_duplicatesRemoved_RGfixed.bam
  do 

  echo " BEGIN ASEREADCOUNT FOR: $SAMPLE"
  bsub -J ASEreadcount -P acc_ISDS -q premium -n 1 -R rusage[mem=16000] -M 80000 -R span[hosts=1] -W 10:00 -o logs/ASEreadcount_%J.stdout -e logs/ASEreadcount_%J.stderr ./3a_ASEReadCounter.sh ${SAMPLE}

done

## EOF ##
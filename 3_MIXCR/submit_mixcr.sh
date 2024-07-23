#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch

for FASTQ_R1_FULL in ${INPUT_DIR}/*R1_filt_001.fastq.gz
  do 

  echo " BEGIN MIXCR ALIGNMENT FOR: $FASTQ_R1_FULL"
  bsub -J MiXCR -P acc_ISDS -q premium -n 8 -R rusage[mem=6000] -R span[hosts=1] -M 384000 -W 96:00 -o logs/MIXCR_%J.stdout -e logs/MIXCR_%J.stderr ./mixcr.sh ${FASTQ_R1_FULL}

done

## EOF ##

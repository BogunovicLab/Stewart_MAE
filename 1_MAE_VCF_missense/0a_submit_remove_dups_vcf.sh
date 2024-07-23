#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_missense

## SAMPLE=${INPUT_DIR}/D108.hard-filtered-noChr.vcf

for SAMPLE in ${INPUT_DIR}/D108_NonSynonymousVariants_PredictedDeleterious-noChr.vcf ${INPUT_DIR}/D113_NonSynonymousVariants_PredictedDeleterious-noChr.vcf ${INPUT_DIR}/D114_NonSynonymousVariants_PredictedDeleterious-noChr.vcf ${INPUT_DIR}/D96F_NonSynonymousVariants_PredictedDeleterious-noChr.vcf
  do 

  echo " BEGIN REMOVE DUPS FOR: $SAMPLE"
  bsub -J remove_dup -P acc_ISDS -q gpu -n 1 -R rusage[mem=16000] -R span[ptile=1] -M 80000 -W 10:00 -o logs/remove-dup_%J.stdout -e logs/remove-dup_%J.stderr ./0a_remove_dups_vcf.sh ${SAMPLE}

done

## EOF ##
#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED

## SAMPLE=${INPUT_DIR}/DON02MAE_F_noDup.vcf.gz

for SAMPLE in ${INPUT_DIR}/D108.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz ${INPUT_DIR}/D113.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz ${INPUT_DIR}/D114.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz ${INPUT_DIR}/D96F.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
  do 

  echo " BEGIN PSEUDOREF FOR: $SAMPLE"
  bsub -J pseudoref -P acc_ISDS -q gpu -n 2 -R rusage[mem=32000] -R span[ptile=2] -M 840000 -W 10:00 -o logs/pseudoref_%J.stdout -e logs/pseudoref_%J.stderr ./0b_create_pseudoRef_hetAllelic.sh ${SAMPLE}

done

## EOF ##
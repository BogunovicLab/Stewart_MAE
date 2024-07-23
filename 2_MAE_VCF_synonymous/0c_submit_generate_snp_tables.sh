#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/HET_ALLELIC

## SAMPLE=${INPUT_DIR}/Het_Allelic_D108.exons.vcf

for SAMPLE in ${INPUT_DIR}/Het_Allelic_D108_noDup.exons.vcf ${INPUT_DIR}/Het_Allelic_D113_noDup.exons.vcf ${INPUT_DIR}/Het_Allelic_D114_noDup.exons.vcf ${INPUT_DIR}/Het_Allelic_D96F_noDup.exons.vcf
  do 

  echo " BEGIN SNP_TABLE FOR: $SAMPLE"
  bsub -J snp_table -P acc_ISDS -q gpu -n 1 -R rusage[mem=16000] -R span[ptile=1] -M 80000 -W 10:00 -o logs/snptable_%J.stdout -e logs/snptable_%J.stderr ./0d_generate_snp_tables.sh ${SAMPLE}

done

## EOF ##
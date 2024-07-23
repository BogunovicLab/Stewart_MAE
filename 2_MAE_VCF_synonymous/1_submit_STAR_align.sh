#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch

for FASTQ_R1_FULL in ${INPUT_DIR}/*R1_filt_001.fastq.gz
  do 

  echo " BEGIN STAR ALIGNMENT FOR: $FASTQ_R1_FULL"
  bsub -J STAR -P acc_ISDS -q premium -n 8 -R rusage[mem=6000] -R span[hosts=1] -M 384000 -W 96:00 -o logs/STAR_%J.stdout -e logs/STAR_%J.stderr ./1_STAR_align.sh ${FASTQ_R1_FULL}

done

## EOF ##

## unzip 
#gunzip /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED/D96F.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
#gunzip /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED/D108.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
#gunzip /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED/D113.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
#gunzip /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED/D114.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz

#bgzip -c D108.HeterozygousSynonymousVariants-noChr_noDup.vcf > D108.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
#bgzip -c D113.HeterozygousSynonymousVariants-noChr_noDup.vcf > D113.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
#bgzip -c D114.HeterozygousSynonymousVariants-noChr_noDup.vcf > D114.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz
#bgzip -c D96F.HeterozygousSynonymousVariants-noChr_noDup.vcf > D96F.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz

# bash 1_submit_STAR_align.sh
#FASTQ_R1_FULL=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch/D41-18_R1_filt_001.fastq.gz

#D108 D113 D41 D68 D96F

## finished very fast
# D108-1
# D113-32
# D113-36
# D113-44
# D41-29

#!/bin/bash

WORK_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/HET_ALLELIC/remove_dup
OUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/SNP_TABLES

module load python
module load gatk
module load vcftools
module load samtools

## SAMPLE_FILE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs/HET_ALLELIC/Het_Allelic_D108.exons.vcf
SAMPLE_FILE_PATH=${1}
SAMPLE_FILE=${SAMPLE_FILE_PATH##*/}
SAMPLE_NAME=${SAMPLE_FILE%%.*}

grep "^#CHROM" ${WORK_DIR}/${SAMPLE_NAME}.exons.vcf | cut -f1-5 > ${OUT_DIR}/${SAMPLE_NAME}.snp_table.txt
grep -v "^#" ${WORK_DIR}/${SAMPLE_NAME}.exons.vcf | sort -V | cut -f1-5 >> ${OUT_DIR}/${SAMPLE_NAME}.snp_table.txt

mv Het_Allelic_D108_noDup.snp_table.txt Het_Allelic_D108.snp_table.txt 
mv Het_Allelic_D113_noDup.snp_table.txt Het_Allelic_D113.snp_table.txt
mv Het_Allelic_D114_noDup.snp_table.txt Het_Allelic_D114.snp_table.txt
mv Het_Allelic_D96F_noDup.snp_table.txt Het_Allelic_D96F.snp_table.txt
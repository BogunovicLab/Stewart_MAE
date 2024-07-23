#!/bin/bash

module load R
module load bedtools

SCRIPTS_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/scripts
REF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE
ASE_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/ASEReadCounter_synonymous-new_batch
SNP_TABLE_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/SNP_TABLES
OUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/COUNTS_synonymous-new_batch

## EX: SAMPLE_FILE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/ASEReadCounter-new_batch/D41-27_ASE.output.table

SAMPLE_FILE_PATH=${1}
SAMPLE_FILE=${SAMPLE_FILE_PATH##*/}
SAMPLE_NAME=${SAMPLE_FILE%%_*}

DONOR=${SAMPLE_NAME%%-*}
VCF_ID=Het_Allelic_D114

Rscript --vanilla ${SCRIPTS_DIR}/counts_to_snps.R -d ${ASE_DIR} -n ${SAMPLE_NAME} -s _ASE.output.table -r ${SAMPLE_NAME} -b ${REF_DIR}/Homo_sapiens.GRCh38.100.EXONS.bed -v ${SNP_TABLE_DIR}/${VCF_ID}.snp_table.txt -o ${OUT_DIR}
#!/bin/bash

module load gatk
module load picard
module load samtools

REF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE
DUP_REMOVE_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_DUP_REMOVED_missense-new_batch
VCF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_missense/HET_ALLELIC/remove_dup
ASE_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/ASEReadCounter_missense-new_batch

## EX: SAMPLE_FILE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_DUP_REMOVED-new_batch/D41-30_duplicatesRemoved_RGfixed.bam

SAMPLE_FILE_PATH=${1}
SAMPLE_FILE=${SAMPLE_FILE_PATH##*/}
SAMPLE_NAME=${SAMPLE_FILE%%_*}
DONOR=${SAMPLE_NAME%%-*}
VCF_ID=Het_Allelic_D114_noDup.exons.vcf.gz

gatk ASEReadCounter --reference ${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I ${DUP_REMOVE_DIR}/${SAMPLE_NAME}_duplicatesRemoved_RGfixed.bam -V ${VCF_DIR}/${VCF_ID} -O ${ASE_DIR}/${SAMPLE_NAME}_ASE.output.table


# Het_Allelic_D96F_noDup.exons.vcf.gz

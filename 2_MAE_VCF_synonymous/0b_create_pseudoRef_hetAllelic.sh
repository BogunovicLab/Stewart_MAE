#!/bin/bash

## directories
REF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE/
VCF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED
HET_VCF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/HET_ALLELIC/
PSEUDOREF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE/PSEUDO_REF/
SCRIPTS=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/scripts

module load python
module load gatk
module load vcftools
module load samtools
module load bcftools

## SAMPLE_FILE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED/D108.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz

SAMPLE_FILE_PATH=${1}
SAMPLE_FILE=${SAMPLE_FILE_PATH##*/}
SAMPLE_NAME=${SAMPLE_FILE%%.*}
VCF_ID=`echo $SAMPLE_NAME`

tabix -p vcf ${VCF_DIR}/${VCF_ID}.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz

python3 ${SCRIPTS}/prepare_reference.py --PSEUDOREF True --HETVCF True \
  --pseudoref_dir ${PSEUDOREF_DIR} \
  --vcf_dir ${HET_VCF_DIR} \
  --ref ${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --name_ind ${VCF_ID} \
  --vcf_ind ${VCF_DIR}/${VCF_ID}.HeterozygousSynonymousVariants-noChr_noDup.vcf.gz \
  --gtf ${REF_DIR}/Homo_sapiens.GRCh38.100.gtf
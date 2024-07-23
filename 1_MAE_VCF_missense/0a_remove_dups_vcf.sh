#!/bin/bash

module load gatk
module load picard
module load samtools
module load bcftools

REF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE
VCF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_missense

# SAMPLE_FILE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_missense/D108_NonSynonymousVariants_PredictedDeleterious-noChr.vcf 

SAMPLE_FILE_PATH=${1}
SAMPLE_FILE=${SAMPLE_FILE_PATH##*/}
SAMPLE_NAME=${SAMPLE_FILE%%_*}
VCF_ID=`echo $SAMPLE_NAME'_NonSynonymousVariants_PredictedDeleterious-noChr'`

## if VCF_ID=DON68MAE_F, cat ${VCF_DIR}/${VCF_ID}.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${VCF_DIR}/${VCF_ID}_PRESORT.vcf
# bgzip -c ${VCF_DIR}/${VCF_ID}_PRESORT.vcf > ${VCF_DIR}/${VCF_ID}.vcf.gz

bgzip -c ${VCF_DIR}/${VCF_ID}.vcf > ${VCF_DIR}/${VCF_ID}.vcf.gz
tabix -p vcf ${VCF_DIR}/${VCF_ID}.vcf.gz
bcftools sort ${VCF_DIR}/${VCF_ID}.vcf.gz -o ${VCF_DIR}/${VCF_ID}_sort.vcf  -Ov

## keep only indels
gatk SelectVariants -R ${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V ${VCF_DIR}/${VCF_ID}_sort.vcf --select-type-to-include SNP -O ${VCF_DIR}/${VCF_ID}_sort_SNPs.vcf
bgzip -c ${VCF_DIR}/${VCF_ID}_sort_SNPs.vcf > ${VCF_DIR}/${VCF_ID}_sort_SNPs.vcf.gz
tabix -p vcf ${VCF_DIR}/${VCF_ID}_sort_SNPs.vcf.gz

## keep normal chrs
bcftools view ${VCF_DIR}/${VCF_ID}_sort_SNPs.vcf.gz --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT -o ${VCF_DIR}/${VCF_ID}_goodChr.vcf

## remove duplicates
bcftools norm --remove-duplicates ${VCF_DIR}/${VCF_ID}_goodChr.vcf -o ${VCF_DIR}/${VCF_ID}_noDup.vcf

#( grep '^#' ${VCF_DIR}/${VCF_ID}_sort.vcf; grep -v "^#" ${VCF_DIR}/${VCF_ID}_sort.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k4,4 | awk -F '\t' 'BEGIN{ prev="";} {key=sprintf("%s\t%s\t%s",$1,$2,$4);if(key==prev) next;print;prev=key;}' ) > ${VCF_DIR}/${VCF_ID}_noDup.vcf

bgzip -c ${VCF_DIR}/${VCF_ID}_noDup.vcf > ${VCF_DIR}/${VCF_ID}_noDup.vcf.gz
gatk IndexFeatureFile -I ${VCF_DIR}/${VCF_ID}_noDup.vcf.gz -O ${VCF_DIR}/${VCF_ID}_noDup.vcf.gz.tbi



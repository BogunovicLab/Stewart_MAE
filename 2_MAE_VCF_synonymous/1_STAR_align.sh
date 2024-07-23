#!/bin/bash

## directories
REF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE/genomeDir/
GTF=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE/Homo_sapiens.GRCh38.100.gtf
FASTQ_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch/
ALIGN_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/STAR_ALIGN_synonymous
VCF_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/VCFs_synonymous/RENAMED/

module load star

FASTQ_R1_FULL=${1}
FASTQ_R1_NAME=${FASTQ_R1_FULL##*/}
SAMPLE_NAME=`echo $FASTQ_R1_NAME | awk -F "_R1"  '{print $1}'`
DONOR=`echo $SAMPLE_NAME | awk -F "-"  '{print $1}'`
FASTQ_R1=${SAMPLE_NAME}_R1_filt_001.fastq.gz
FASTQ_R2=${SAMPLE_NAME}_R2_filt_001.fastq.gz
VCF_ID=D114.HeterozygousSynonymousVariants-noChr_noDup.vcf

mkdir -p ${ALIGN_DIR}/${DONOR}/${SAMPLE_NAME}/

STAR \
    --readFilesIn ${FASTQ_DIR}/${FASTQ_R1} ${FASTQ_DIR}/${FASTQ_R2} \
    --outFileNamePrefix ${ALIGN_DIR}/${DONOR}/${SAMPLE_NAME}/ \
    --runThreadN 8 \
    --outSAMtype BAM SortedByCoordinate  \
    --outSAMattributes vG vA NH HI AS nM NM \
    --genomeDir ${REF_DIR} \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --sjdbGTFfile $GTF \
    --varVCFfile ${VCF_DIR}/${VCF_ID} \
    --waspOutputMode SAMtag \
    --twopassMode Basic  \
    --quantMode GeneCounts

echo "FINISHED STAR ALIGNMENT FOR ${SAMPLE_NAME}"
echo "VCF: ${VCF_ID}"

#cp D96F.hard-filtered-noChr_noDup.vcf RENAMED/
#cp D114.hard-filtered-noChr_noDup.vcf RENAMED/
#cp D113.hard-filtered-noChr_noDup.vcf RENAMED/
#cp D108.hard-filtered-noChr_noDup.vcf RENAMED/
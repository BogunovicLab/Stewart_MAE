#!/bin/bash

module load samtools
module load gatk
module load picard

BAM_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_ALIGNED_missense-new_batch
WASP_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_WASP_FILT_missense-new_batch
DUP_REMOVE_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_DUP_REMOVED_missense-new_batch

## EX: SAMPLE_FILE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_ALIGNED-new_batch/D108-10_Aligned.sortedByCoord.out.bam

SAMPLE_FILE_PATH=${1}
SAMPLE_FILE=${SAMPLE_FILE_PATH##*/}
SAMPLE_NAME=${SAMPLE_FILE%%_*}

## only keep reads that passed WASP filtering
(samtools view -H ${BAM_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam; samtools view ${BAM_DIR}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam | grep --color 'vW:i:1') | samtools view -bS - > ${WASP_DIR}/${SAMPLE_NAME}_WASPfilt.bam

## remove duplicates
java -jar $PICARD MarkDuplicates I=${WASP_DIR}/${SAMPLE_NAME}_WASPfilt.bam O=${DUP_REMOVE_DIR}/${SAMPLE_NAME}_duplicatesRemoved.bam M=${DUP_REMOVE_DIR}/${SAMPLE_NAME}_duplicatesRemoved.marked_dup_metrics.txt REMOVE_DUPLICATES=true

## fix read groups so downstream programs don't complain
java -jar $PICARD AddOrReplaceReadGroups I=${DUP_REMOVE_DIR}/${SAMPLE_NAME}_duplicatesRemoved.bam O=${DUP_REMOVE_DIR}/${SAMPLE_NAME}_duplicatesRemoved_RGfixed.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20


# samtools view ${WASP_DIR}/D02c-Clone-E10a4_WASPfilt.bam | grep --color 'vW:i:1'

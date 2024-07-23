#!/bin/bash

fastqDir=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch
inDir=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/STAR_ALIGN_missense
outDir=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_ALIGNED_missense-new_batch

for FASTQ_R1_FULL in ${fastqDir}/*R1_filt_001.fastq.gz
	do

	## EX: FASTQ_R1_FULL=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch/D96F-16_R1_filt_001.fastq.gz

	FASTQ_R1_NAME=${FASTQ_R1_FULL##*/}
	SAMPLE_NAME=`echo $FASTQ_R1_NAME | awk -F "_R1"  '{print $1}'`
	DONOR=`echo $SAMPLE_NAME | awk -F "-"  '{print $1}'`

	cp ${inDir}/${DONOR}/${SAMPLE_NAME}/Aligned.sortedByCoord.out.bam ${inDir}/${DONOR}/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam
	mv ${inDir}/${DONOR}/${SAMPLE_NAME}/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam ${outDir}

done

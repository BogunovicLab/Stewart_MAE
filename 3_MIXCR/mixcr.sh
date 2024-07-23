#!/bin/bash

ml mixcr
## MiXCR v3.0.13

INPUT_DIR="/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch/"
OUTPUT_DIR="/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/MIXCR/"

# parameters
species="hs"
starting_material="rna"
threads="24"

#SAMPLE_PATH=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/TRIMMED_FASTQS-new_batch/D108-10_R1_filt_001.fastq.gz
SAMPLE_PATH=${1}
SAMPLE_NAME=$(basename ${SAMPLE_PATH} | sed 's/_R1_filt_001.fastq.gz//')

# Run mixcr analyze shotgun
mixcr analyze shotgun \
        -s ${species} \
        -t ${threads} \
        --starting-material ${starting_material} \
        ${SAMPLE_PATH} \
        ${SAMPLE_PATH/_R1_/_R2_} \
        ${OUTPUT_DIR}/${SAMPLE_NAME}

echo "finished MiXCR alignment for: ${SAMPLE_NAME}"
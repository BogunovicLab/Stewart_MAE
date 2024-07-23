#!/bin/bash

INPUT_DIR=/sc/arion/projects/ISDS/haleyr/RMAE_OJAY/DATA/BAMS_ALIGNED_missense-new_batch

for SAMPLE in ${INPUT_DIR}/*_Aligned.sortedByCoord.out.bam
  do 

  echo " BEGIN WASP, DEDUP FOR: $SAMPLE"
  bsub -J WASP -P acc_ISDS -q premium -n 6 -R rusage[mem=10000] -R span[hosts=1] -M 320000 -W 10:00 -o logs/WASP_%J.stdout -e logs/WASP_%J.stderr ./2_WASP_and_removeDup_filters.sh ${SAMPLE}

done

## EOF ##
#!/bin/bash

awk '$3=="exon" && ($1 ~ /^[1-9]/ || $1 == "X" || $1 == "Y")' /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE/Homo_sapiens.GRCh38.100.gtf | cut -f1,4,5,9 | awk -F $'\t' 'BEGIN {OFS = FS} {split($4, a, "gene_id \""); split(a[2], b, "\""); print $1, $2-1, $3, b[1]}' > /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/REFERENCE/Homo_sapiens.GRCh38.100.EXONS.bed

#!/bin/bash

# run qualimap-bamqc
/opt/PepPrograms/qualimap_v2.2.1/qualimap bamqc -bam $1.sort.bam -outdir QC/$1_qualimap/bamqc/

# run qualimap-rnaseq
/opt/PepPrograms/qualimap_v2.2.1/qualimap rnaseq -bam $1.sort.bam -outdir QC/$1_qualimap/rnaseq/ -pe -gtf $2

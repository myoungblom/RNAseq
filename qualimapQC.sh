#!/bin/bash

# run qualimap-bamqc
/opt/PepPrograms/qualimap_v2.2.1/qualimap bamqc -bam $1.sort.bam -outdir $1_bamqc --java-mem-size=2G

# run qualimap-rnaseq
/opt/PepPrograms/qualimap_v2.2.1/qualimap rnaseq -bam $1.sort.bam -outdir $1_rnaseq -pe -gtf $2 --java-mem-size=2G

tar -zcvf $1.qualimap.tar.gz $1_bamqc/ $1_rnaseq/

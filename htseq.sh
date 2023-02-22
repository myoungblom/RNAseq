#!/bin/bash

python3 -m HTSeq.scripts.count -f bam -r pos -t gene -i gene_name --nonunique none -s reverse $1.sort.bam $2

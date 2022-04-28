#!/bin/bash

python3 -m HTSeq.scripts.count -f bam -r pos -t gene -i gene_name --nonunique all -s no $1.sort.bam $2

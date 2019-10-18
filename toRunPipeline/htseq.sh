#!/bin/bash

# run HTseq Counts

python3 -m HTSeq.scripts.count -f bam -r pos -t gene -i Name --nonunique none -s reverse $1.sort.bam $2

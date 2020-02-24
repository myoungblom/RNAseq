#!/bin/bash

python3 -m HTSeq.scripts.count -f bam -r pos -t gene -i Name --nonunique none -s no $1.sort.bam $2

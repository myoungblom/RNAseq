#!/bin/bash

# check that counts file exists
file1='find $1 -type f -size +$2 -print'

if [ -z "$file1" ]; then
	exit 1
fi

# create reduced counts file
file2="${1/.txt/_reduced.txt}"
cut -f 1,3 "$1" | head -n -5 > "HTSeqCounts/$file2"

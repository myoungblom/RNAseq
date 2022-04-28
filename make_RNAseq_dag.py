#!/usr/bin/env python

#####
# This script takes an input file containing metadata from RNA sequencing and
# will output a dag file for each sample.
#####

# Format for input file: one sample name per line

import sys
import os
import argparse
from string import Template

def get_args():
    parser = argparse.ArgumentParser(description='Pipeline for RNAseq data preparation, alignment and processing')
    parser.add_argument('input', help='File describing read data information')
    parser.add_argument('dagtemplate', help='Dag template')
    parser.add_argument('reference', help='Path to reference fasta file')
    parser.add_argument('gff', help='Path to GFF annotation file')
    parser.add_argument('gtf', help='Path to GTF annotation file')
    
    return parser.parse_args()

args = get_args()

with open(args.dagtemplate, 'r') as template_file:
    template = Template(template_file.read())

with open(args.input, 'r') as infile:
    outFile = args.input.split('.')[0]+'_topLevel.dag'
    with open(outFile, 'w') as output:
        for line in infile:
            line = line.strip('\n')
            inputList = line.split('\t')
            variableMap = {}
            variableMap['ref'] = args.reference
            variableMap['annot_gff'] = args.gff
            variableMap['annot_gtf'] = args.gtf
            variableMap['run'] = inputList[0]
            with open('{0}_RNAseq.dag'.format(variableMap['run']), 'w') as dagfile:
                out = template.substitute(variableMap)
                dagfile.write(out)

            output.write('SPLICE {0} {0}_RNAseq.dag\n'.format(variableMap['run']))


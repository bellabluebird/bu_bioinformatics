#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help='a GFF file', dest="input", required=True)
parser.add_argument("-o", "--output", help='output file with region', dest="output", required=True)

args = parser.parse_args()

import re
id_2_name = {}

genename = r'gene_name\s([^;]*)'
geneid = r'gene_id\s([^;]*)'

with open(args.input, 'r') as r:
    for line in r:
        if line.startswith('#'):
            continue
        gene_name = re.search(genename,line)
        gene_id = re.search(geneid,line)

        if gene_id.group().split('"')[1] in id_2_name:
            continue
        else: 
            id_2_name[gene_id.group().split('"')[1]] = gene_name.group().split('"')[1]

with open(args.output, 'wt') as w: 
    for k,v, in id_2_name.items():
        w.write('{}\t{}\n'.format(k,v))
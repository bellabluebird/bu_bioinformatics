#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", required=True)
parser.add_argument("-o", "--output", dest="output", required=True)
parser.add_argument("-s", "--subset", dest="subset", required=True)

args = parser.parse_args()

with open(args.input, 'rt') as handle, open(args.output, 'wt') as w:
    for rec in SeqIO.parse(handle, 'fasta'):
        if rec.id == args.subset:
            SeqIO.write(rec, w, "fasta")
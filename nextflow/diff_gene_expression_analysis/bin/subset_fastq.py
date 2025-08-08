#!/usr/bin/env python

import argparse
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help='fastq file',dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file name and path',dest="output", required=True)
parser.add_argument("-s", "--srr", help='SRR ID',dest="srr", required=True)
parser.add_argument("-n", "--new", help='New name',dest="new", required=True)

args = parser.parse_args()

def update_element(slot, old, new):
    return(slot.replace(old, new))

out = []
with gzip.open(args.input, 'rt') as handle, open(args.output, "w") as output_handle:
    for rec in SeqIO.parse(handle, 'fastq'):
        repl = {'id': update_element(rec.id, args.srr, args.new), 'name': update_element(rec.name, args.srr, args.new), 'description': update_element(rec.description, args.srr, args.new)}
        vars(rec).update(repl)
        SeqIO.write(rec, output_handle, "fastq")
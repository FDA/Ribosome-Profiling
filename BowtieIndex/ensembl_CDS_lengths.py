#!/usr/bin/env python

from Bio import SeqIO
import argparse as ap

parser = ap.ArgumentParser(description = "Parse CDS lengths")
parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
args = parser.parse_args()

cds_lengths = {}

infh = open("%s_all_CDS.fa" % args.dataset, "r")

for record in SeqIO.parse(infh, "fasta"):
    cds_lengths[record.id] = len(record.seq)
    

outfh = open("%s_ensembl_CDS_lengths.tsv" % args.dataset, "w")

for key in cds_lengths:
    outfh.write(key + "\t" + str(cds_lengths[key]) + "\n")

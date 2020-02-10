#!/usr/bin/env python

import re
import argparse as ap 

parser = ap.ArgumentParser(description = "Parse transcript lengths")
parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
args = parser.parse_args()


pattern = re.compile(r'(?=transcript_id=([^;]*))')

lengths = {}

infh = open("%s_all_annotations.gff3" % args.dataset, "r")
for line in infh:
    try:
        m = pattern.search(line)
        lengths[m.group(1)] = 0
    except AttributeError as e:
        continue
infh.close()

infh = open("%s_all_annotations.gff3" % args.dataset, "r")

for line in infh:
    s = line.split("\t")
    #print s
    
    if len(s) > 5:
        type = s[2]
        if type == "exon":
            try:
                m = pattern.search(line)
                tid = m.group(1)
                start, end = int(s[3]), int(s[4])
                lengths[tid] += (end-start)+1
                #print start, end
            except AttributeError as e:
                print line
                
                
outfh = open("%s_ensembl_transcript_lengths.tsv" % args.dataset, "w")
for key in lengths.keys():
    outfh.write(key + "\t" + str(lengths[key]) + "\n")
outfh.close()
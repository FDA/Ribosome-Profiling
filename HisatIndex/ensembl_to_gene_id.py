#!/usr/bin/env python

import re
import argparse as ap 

parser = ap.ArgumentParser(description = "Parse gene names")
parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
args = parser.parse_args()

infh = open("%s_all_annotations.gff3" % args.dataset, "r")


p1 = re.compile(r'(?=transcript_id=([^;]*))')
p2 = re.compile(r'(?=gene_name=([^;]*))')

convert = {}

for line in infh:
    if line[0] == "#": 
        continue
    m1 = p1.search(line)
    m2 = p2.search(line)
    
    try:
        tid = m1.group(1)
        gid = m2.group(1)
        convert[tid] = gid
        
    except:
        print(line)
        
        
outfh = open("%s_ensembl_to_gene_id.tsv" % args.dataset, "w")
for key in convert.keys():
    outfh.write(key + "\t" + str(convert[key]) + "\n")
outfh.close()

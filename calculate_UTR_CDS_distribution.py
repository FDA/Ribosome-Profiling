#!/usr/bin/env python

import pysam
import sys
import os
import time
import collections
import argparse
start_time = time.time()


TypeDict = collections.OrderedDict([("5UTR", 0), ("3UTR", 0), ("CDS", 0)])

parser = argparse.ArgumentParser(description = "Calculate the reads aligning to each codon of a gene of interest.")
parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
args = parser.parse_args()

type_data = {}

data_files = sorted([x for x in os.listdir("./Annotated_reads/%s/" % args.dataset) if x[-4:] == ".bam"])
base_names = []

for f in data_files:
    base_name = "_".join(f.split(".")[0].split("_")[1:])  # key for all_reads
    base_names.append(base_name)
    type_data[base_name] = {}
    
    for i in range(20, 126):
        type_data[base_name][i] = TypeDict.copy()

    read_totals = {}
    base_name = "_".join(f.split(".")[0].split("_")[1:])
    tot_infh = open("./Annotated_totals_per_read/%s/%s_read_counts.tsv" % (args.dataset, base_name), "r")
    tot_infh.next()  # skip header
    for line in tot_infh:
        s = line.rstrip().split("\t")
        read_totals[s[0]] = float(s[1])
    tot_infh.close()
    
    bamfile = pysam.AlignmentFile("./Annotated_reads/%s/%s" % (args.dataset, f), "rb")
    for read in bamfile.fetch(until_eof = True):
        t = read.get_tag("ut")
        length = read.query_length
        read_nh = read_totals[read.qname]
        type_data[base_name][length][t] += (1/read_nh)
        
        

if not os.path.exists("./output_UTR_dist_by_length"):
    os.mkdir("./output_UTR_dist_by_length")
if not os.path.exists("./output_UTR_dist_by_length/%s" % args.dataset):
    os.mkdir("./output_UTR_dist_by_length/%s" % args.dataset)
    
    
for key in sorted(type_data.keys()):    
    outfh = open("./output_UTR_dist_by_length/%s/%s_%s_UTR_CDS_distribution.csv" % (args.dataset, args.dataset, key), "w")
    header = ",".join(["Length", "CDS", "5UTR", "3UTR", "CDS %", "5UTR %", "3UTR %"]) + "\n"
    outfh.write(header)
    
    all_line = ["All"]
    tot_cds, tot_fputr, tot_tputr = 0, 0, 0
    for length in sorted(type_data[key].keys()):
        tot_cds += type_data[key][length]["CDS"]
        tot_fputr += type_data[key][length]["5UTR"]
        tot_tputr += type_data[key][length]["3UTR"]
        
    all_line.append(str(tot_cds))
    all_line.append(str(tot_fputr))
    all_line.append(str(tot_tputr))
    all_line.append(str(round((float(100*tot_cds)/(tot_cds+tot_fputr+tot_tputr)),2)))
    all_line.append(str(round((float(100*tot_fputr)/(tot_cds+tot_fputr+tot_tputr)),2)))
    all_line.append(str(round((float(100*tot_tputr)/(tot_cds+tot_fputr+tot_tputr)),2)))
    
    outfh.write(",".join(all_line) + "\n")
    
    for length in sorted(type_data[key].keys()):
        cds = type_data[key][length]["CDS"]  # cds count
        fputr = type_data[key][length]["5UTR"]  # 5' utr count
        tputr = type_data[key][length]["3UTR"]  # 3' utr count
        if sum([cds, fputr, tputr]) == 0:
            pcds, pfputr, ptputr = 0, 0, 0
        else:
            pcds = round((float(100*cds)/(cds+fputr+tputr)),2)
            pfputr = round((float(100*fputr)/(cds+fputr+tputr)),2)  # % 5' utr
            ptputr = round((float(100*tputr)/(cds+fputr+tputr)),2)  # % 3' utr
        outline = [str(length), str(cds), str(fputr), str(tputr), str(pcds), str(pfputr), str(ptputr)]
        outfh.write(",".join(outline) + "\n")
        
        
    outfh.close()
    
    

typelist = ["5UTR", "% 5UTR", "CDS", "% CDS", "3UTR", "% 3UTR"]
outfh = open("./output_UTR_dist_by_length/%s/%s_all_UTR_CDS_distribution.csv" % (args.dataset, args.dataset), "w")
header = ["Sample", "Type"]
for i in range(20, 33):
    header.append(str(i))
outfh.write(",".join(header) + "\n")

for key in sorted(type_data.keys()):
    for dtype in typelist:
        outlist = [key, dtype]
        for length in sorted(type_data[key].keys()):
            if dtype == "5UTR":
                outlist.append(str(type_data[key][length]["5UTR"]))
            
            elif dtype == "% 5UTR":
                try:
                    val = round((float(100*type_data[key][length]["5UTR"])/sum(type_data[key][length].values())),2)
                except ZeroDivisionError as e:
                    val = 0
                outlist.append(str(val))
            
            elif dtype == "CDS":
                outlist.append(str(type_data[key][length]["CDS"]))
                
            elif dtype == "% CDS":
                try:
                    val = round((float(100*type_data[key][length]["CDS"])/sum(type_data[key][length].values())),2)
                except ZeroDivisionError as e:
                    val = 0
                outlist.append(str(val))
                
            elif dtype == "3UTR":
                outlist.append(str(type_data[key][length]["3UTR"]))
                
            elif dtype == "% 3UTR":
                try:
                    val = round((float(100*type_data[key][length]["3UTR"])/sum(type_data[key][length].values())),2)
                except ZeroDivisionError as e:
                    val = 0
                outlist.append(str(val))
                
        outfh.write(",".join(outlist) + "\n")
outfh.close()

ribo_base_names = []
for base_name in base_names:
    if base_name[1] == "R":
        ribo_base_names.append(base_name)


outfh = open("./output_UTR_dist_by_length/%s/%s_compiled_ribo_UTR_distribution.csv" % (args.dataset, args.dataset), "w")
header1 = ["Length:"]
for i in range(20,36):
    header1.append(str(i))
    for j in range(0,len(ribo_base_names)-1):
        header1.append("")
outfh.write(",".join(header1) + "\n")

header2 = ["Type"]
for i in range(20,36):
    for ribo_base_name in ribo_base_names:
        header2.append(ribo_base_name)
outfh.write(",".join(header2) + "\n")

fivePrimeUTRline = ["% 5UTR"]
for i in range(20,36):
    for ribo_base_name in ribo_base_names:
        try:
            fivePrimeUTRline.append(str(round(float(100*type_data[ribo_base_name][i]["5UTR"])/sum(type_data[ribo_base_name][i].values()),2)))
        except ZeroDivisionError as e:
            fivePrimeUTRline.append("0")
outfh.write(",".join(fivePrimeUTRline) + "\n")

CDSline = ["% CDS"]
for i in range(20,36):
    for ribo_base_name in ribo_base_names:
        try:
            CDSline.append(str(round(float(100*type_data[ribo_base_name][i]["CDS"])/sum(type_data[ribo_base_name][i].values()),2)))
        except ZeroDivisionError as e:
            CDSline.append("0")
outfh.write(",".join(CDSline) + "\n")

threePrimeUTRLine = ["% 3UTR"]
for i in range(20,36):
    for ribo_base_name in ribo_base_names:
        try:
            threePrimeUTRLine.append(str(round(float(100*type_data[ribo_base_name][i]["3UTR"])/sum(type_data[ribo_base_name][i].values()),2)))
        except ZeroDivisionError as e:
            threePrimeUTRLine.append("0")
outfh.write(",".join(threePrimeUTRLine) + "\n")

outfh.close()


    
    

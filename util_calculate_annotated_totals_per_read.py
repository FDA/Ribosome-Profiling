#!/usr/bin/env python

import sys
import os
import math
import subprocess
import argparse
import time
from collections import OrderedDict
import pysam


def create_output_folders(dataset):
    if not os.path.exists("./Annotated_totals_per_read"):
        os.mkdir("./Annotated_totals_per_read")
    if not os.path.exists("./Annotated_totals_per_read/%s" % dataset):
        os.mkdir("./Annotated_totals_per_read/%s" % dataset)

def calculate_Annotated_totals_per_read(dataset, base_name):
    
    read_totals = {}
    
    bamfile = pysam.AlignmentFile("./Annotated_reads/%s/annotated_%s.bam" % (dataset, base_name), "rb")
    for read in bamfile.fetch():
        # just populating dict the first time
        read_totals[read.query_name] = 0
    bamfile.close()
    
    
    bamfile = pysam.AlignmentFile("./Annotated_reads/%s/annotated_%s.bam" % (dataset, base_name), "rb")
    for read in bamfile.fetch():
        read_totals[read.query_name] += 1
    bamfile.close()
    
    outfh = open("./Annotated_totals_per_read/%s/%s_read_counts.tsv" % (dataset, base_name), "w")
    outfh.write("Read\tCount\n")
    for name in read_totals.keys():
        outfh.write(name + "\t" + str(read_totals[name]) + "\n")
    outfh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate the number of alignments each read has after filtering with.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Annotated_reads folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    
    data_files = sorted(os.listdir("./Annotated_reads/%s/" % args.dataset))  
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[1:])
        calculate_Annotated_totals_per_read(args.dataset, base_name)


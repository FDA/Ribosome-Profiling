#!/usr/bin/env python

import sys
import os
import re
import subprocess
import argparse
import time
from collections import OrderedDict
import pysam

OffsetTally = OrderedDict()
for i in range(-30,21):
    OffsetTally[i] = 0
    
    
def create_output_folders(dataset):
    if not os.path.exists("./output_stop_codon_p_site_offset_by_length"):
	    os.mkdir("./output_stop_codon_p_site_offset_by_length")
    if not os.path.exists("./output_stop_codon_p_site_offset_by_length/%s" % args.dataset):
	    os.mkdir("./output_stop_codon_p_site_offset_by_length/%s" % args.dataset)
    
    
def get_start_end_coords(dataset):
    coordinates_dict = {}  # key: transcript id, value: dict w/ keys "start" and "end"
    infh = open("%s_transcript_CDS_start_end_locations.txt" % dataset, "r")
    infh.next()  # skip header line
    for line in infh:
        s = line.split("\t")
        coordinates_dict[s[0]] = {"start":int(s[1])-1, "end":int(s[2])}  # subtract 1 from start bc the values are 1-indexed, end is inclusive of last base so it balances the -1 shift and python's exclusivity already
        
    #for tid in coordinates_dict.keys():
        #print tid, coordinates_dict[tid]["start"], coordinates_dict[tid]["end"]
        
    return coordinates_dict

def process_sample_p_site_offset(dataset, base_name, coordinates_dict):
    offset_data = OrderedDict()
    for i in range(18,36):
        offset_data[i] = OffsetTally.copy()
     
    
    bamfile = pysam.AlignmentFile("./Tophat_output/%s/%s/%s_sorted.bam" % (dataset, base_name, base_name), "rb")
    
    d = {}
    
    for read in bamfile.fetch():
        if read.query_length > 35:
            continue 
        offset = read.reference_start - (coordinates_dict[read.reference_name]["end"]-3)  # negative offset = upstream
        #offset = read.reference_start - coordinates_dict[read.reference_name]["start"]  # negative offset = upstream
        
        if (offset < -30) or (offset > 20):
            continue
            
        if abs(offset) > read.query_length-3:
            continue

        read_nh = int(read.get_tag("NH"))
        offset_data[read.query_length][offset] += float(1)/read_nh
    
    bamfile.close()
    #print d
    
    outfh = open("./output_stop_codon_p_site_offset_by_length/%s/%s_stop_codon_p_site_offsets.csv" % (dataset, base_name), "w")
    header = "Read Length"
    for l in OffsetTally.keys():
        header += "," + str(l)
    outfh.write(header + "\n")
    
    for i in offset_data.keys():
        outstr = str(i)
        for j in OffsetTally.keys():
            outstr += "," + str(offset_data[i][j])
        outstr += "\n"
        outfh.write(outstr)
        
    outfh.close()
    
    print "Finished with %s at " % base_name, time.ctime()
     
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate the p-site offset relative to the stop codon  of the reads in a dataset after adapter trimming, rRNA+tRNA decontamination, and alignment.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the TophatOutput folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    coordinates_dict = get_start_end_coords(args.dataset)
    
    data_folders = sorted(os.listdir("./Tophat_output/%s/" % args.dataset))
    
    for f in data_folders:#[0:1]:
        base_name = f
        process_sample_p_site_offset(args.dataset, base_name, coordinates_dict)

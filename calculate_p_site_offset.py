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
for i in range(-20,21):
    OffsetTally[i] = 0
    
    
def create_output_folders(dataset):
    if not os.path.exists("./output_p_site_offset_by_length"):
	    os.mkdir("./output_p_site_offset_by_length")
    if not os.path.exists("./output_p_site_offset_by_length/%s" % args.dataset):
	    os.mkdir("./output_p_site_offset_by_length/%s" % args.dataset)
    
    
def get_start_end_coords(dataset):
    #For coordinates_dict... key: transcript id, value: dict w/ keys "start" and "end"
    coordinates_dict = {}
    infh = open("%s_transcript_CDS_start_end_locations.txt" % dataset, "r")
    # skip header line
    next(infh)
    for line in infh:
        s = line.split("\t")
        # subtract 1 from start bc the values are 1-indexed, end is inclusive of last base so it balances the -1 shift
        coordinates_dict[s[0]] = {"start":int(s[1])-1, "end":int(s[2])}
                
    return coordinates_dict

def process_sample_p_site_offset(dataset, base_name, coordinates_dict):
    offset_data = OrderedDict()
    for i in range(18,36):
        offset_data[i] = OffsetTally.copy()
    
    bamfile = pysam.AlignmentFile("./Hisat_output/%s/%s/%s_sorted.bam" % (dataset, base_name, base_name), "rb")
    
    for read in bamfile.fetch():
        if read.query_length > 35:
            continue 
        offset = read.reference_start - coordinates_dict[read.reference_name]["start"]
        
        if abs(offset) > 20:
            continue
            
        
        read_nh = int(read.get_tag("NH"))
        offset_data[read.query_length][offset] += float(1)/read_nh
    
    bamfile.close()
    
    outfh = open("./output_p_site_offset_by_length/%s/%s_p_site_offsets.csv" % (dataset, base_name), "w")
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
    
    print("Finished with %s at " %base_name, time.ctime())
    
    return offset_data
    
    
def write_combined_file(dataset, base_names, offset_data_container):
    ribo_base_names = []
    for base_name in base_names:
        if base_name[1] == "R":
            ribo_base_names.append(base_name)
    
    
    outfh = open("./output_p_site_offset_by_length/%s/%s_combined_ribo_p_site_offsets.csv" % (dataset, dataset), "w")

    header1 = ["Read Length:"]
    for i in range(20, 36):
        header1.append(str(i))
        for j in range(0, len(ribo_base_names)-1):
            header1.append("")
    outfh.write(",".join(header1) + "\n")
    
    header2 = ["Offset"]
    for i in range(20,36):
        for ribo_base_name in ribo_base_names:
            header2.append(str(ribo_base_name))
    outfh.write(",".join(header2) + "\n")
    
    
    for i in range(-20,21):
        outline = [str(i)] 
        for j in range(20, 36):
            for ribo_base_name in ribo_base_names:
                outline.append(str(offset_data_container[ribo_base_name][j][i]))
                
        outfh.write(",".join(outline) + "\n")
        
    outfh.close()
    print("Finished with combined offset file at", time.ctime())
            
        
     
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate the length distribution of the reads in a dataset after adapter trimming, rRNA+tRNA decontamination, and alignment.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the HisatOutput folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    coordinates_dict = get_start_end_coords(args.dataset)
    
    offset_data_container = {}
    
    data_folders = sorted(os.listdir("./Hisat_output/%s/" % args.dataset))
    base_names = []
    
    for f in data_folders:
        base_name = f
        base_names.append(base_name)
        offset_data_container[base_name] = process_sample_p_site_offset(args.dataset, base_name, coordinates_dict)
        
    write_combined_file(args.dataset, base_names, offset_data_container)



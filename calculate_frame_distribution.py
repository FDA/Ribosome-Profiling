#!/usr/bin/env python

import pysam
import sys
import os
import time
import collections
import argparse
start_time = time.time()


FrameDict = collections.OrderedDict([(0, 0), (1, 0), (2, 0)])

parser = argparse.ArgumentParser(description = "Calculate the distribution of reading frames in a dataset.")
parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
args = parser.parse_args()

frame_data = {}

data_files = sorted([x for x in os.listdir("./Annotated_reads/%s/" % args.dataset) if x[-4:] == ".bam"])
base_names = []

for f in data_files:
    base_name = "_".join(f.split(".")[0].split("_")[1:])  # key for all_reads
    base_names.append(base_name)
    frame_data[base_name] = {}
    
    for i in range(20, 126):
        frame_data[base_name][i] = FrameDict.copy()


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
        f = read.get_tag("rf")
        if f == -1:
            continue  # skip the UTR reads that don't have a frame
        length = read.query_length
        read_nh = read_totals[read.qname]
        frame_data[base_name][length][f] += float(1)/read_nh
        
        

if not os.path.exists("./output_frame_dist_by_length"):
    os.mkdir("./output_frame_dist_by_length")
if not os.path.exists("./output_frame_dist_by_length/%s" % args.dataset):
    os.mkdir("./output_frame_dist_by_length/%s" % args.dataset)
    
    
for key in sorted(frame_data.keys()):    
    outfh = open("./output_frame_dist_by_length/%s/%s_%s_frame_distribution.csv" % (args.dataset, args.dataset, key), "w")
    header = ",".join(["Length", "+0", "+1", "+2", "+0 %", "+1 %", "+2 %"]) + "\n"
    outfh.write(header)
    
    all_line = ["All"]
    tot_frame_zero, tot_frame_one, tot_frame_two = 0, 0, 0
    for length in sorted(frame_data[key].keys()):
        tot_frame_zero += frame_data[key][length][0]
        tot_frame_one += frame_data[key][length][1]
        tot_frame_two += frame_data[key][length][2]
        
    all_line.append(str(tot_frame_zero))
    all_line.append(str(tot_frame_one))
    all_line.append(str(tot_frame_two))
    all_line.append(str(round((float(100*tot_frame_zero)/(tot_frame_zero+tot_frame_one+tot_frame_two)),2)))
    all_line.append(str(round((float(100*tot_frame_one)/(tot_frame_zero+tot_frame_one+tot_frame_two)),2)))
    all_line.append(str(round((float(100*tot_frame_two)/(tot_frame_zero+tot_frame_one+tot_frame_two)),2)))
    
    outfh.write(",".join(all_line) + "\n")
    
    for length in sorted(frame_data[key].keys()):
        frame_zero = frame_data[key][length][0]
        frame_one = frame_data[key][length][1]
        frame_two = frame_data[key][length][2]
        if sum([frame_zero, frame_one, frame_two]) == 0:
            pfzero, pfone, pftwo = 0, 0, 0
        else:
            pfzero = round((float(100*frame_zero)/(frame_zero+frame_one+frame_two)),2)
            pfone = round((float(100*frame_one)/(frame_zero+frame_one+frame_two)),2)
            pftwo = round((float(100*frame_two)/(frame_zero+frame_one+frame_two)),2)
        outline = [str(length), str(frame_zero), str(frame_one), str(frame_two), str(pfzero), str(pfone), str(pftwo)]
        outfh.write(",".join(outline) + "\n")
    
    outfh.close()
        
    print "Finished with %s frame distribution at " % key, time.ctime()

print "Now starting compiled file at", time.ctime(), "..."    

outfh = open("./output_frame_dist_by_length/%s/%s_frame_distribution_compiled.csv" % (args.dataset, args.dataset), "w")
num_files = len(data_files)
first_header = [""] # space for frame column
for i in range(20,126):
    first_header.append(str(i))  # number...
    for j in range(0,num_files-1):
        first_header.append("")  # ...followed by n-1 spaces until next number
outfh.write(",".join(first_header) + "\n")

second_header = ["Frame"]
for i in range(20,126):
    for f in base_names:
        second_header.append(f)
outfh.write(",".join(second_header) + "\n")

for frame in ["0", "1", "2"]:
    frame_line = [frame] 
    for i in range(20,126):
        for base_name in base_names:
            if (frame_data[base_name][i][0]+frame_data[base_name][i][1]+frame_data[base_name][i][2]) == 0:
                number = 0
            else:
                number = round((float(100*frame_data[base_name][i][int(frame)])/(frame_data[base_name][i][0]+frame_data[base_name][i][1]+frame_data[base_name][i][2])),2)
            frame_line.append(str(number))
            
    outfh.write(",".join(frame_line) + "\n")
outfh.close()

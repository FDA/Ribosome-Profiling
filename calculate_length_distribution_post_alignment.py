#!/usr/bin/env python

from Bio import SeqIO
import argparse
import pysam
import time
import collections
import sys
import os


def calculate_length_distribution(dataset):
    # Create a template container for the tallies of read lengths to be stored #
    LengthTallyDict = collections.OrderedDict()
    for i in range(0, 126):
        LengthTallyDict[i] = 0
    
    # Get list of input files, in .bam format #
    data_folders = sorted(os.listdir("./Tophat_output/%s/" % args.dataset))  # use the accepted_hits.bam folder in each file
    
    # Create folders for the output to be placed #
    if not os.path.exists("./output_length_dist_post_alignment"):
        os.mkdir("./output_length_dist_post_alignment")
    if not os.path.exists("./output_length_dist_post_alignment/%s" % args.dataset):
        os.mkdir("./output_length_dist_post_alignment/%s" % args.dataset)
        
    
    # Create output file (.csv) to place the counts of reads with each length
    # Header looks like Sample,18,19,20,21 ... 
    # First row looks like 1R_CAGATC_L002_R1_001,12345,23456, ... 
    outfh_len = open("./output_length_dist_post_alignment/%s/%s_post_alignment_length_counts.csv" % (args.dataset, args.dataset), "w")
    header_len = "Sample"
    for i in LengthTallyDict.keys():
        header_len += "," + str(i)
    outfh_len.write(header_len + "\n")


    # Create output file (.csv) to place the percents of reads with each length
    # Header looks like Sample,% Reads Len 18,% Reads Len 19, ... 
    # First row looks like 1R_CAGATC_L002_R1_001,12.345,23.456, ... 
    outfh_lenperc = open("./output_length_dist_post_alignment/%s/%s_post_alignment_length_percents.csv" % (args.dataset, args.dataset), "w")
    header_len = "Sample"
    for i in LengthTallyDict.keys():
        header_len += ",% Reads Len " + str(i)
    outfh_lenperc.write(header_len + "\n")
    
    
    for f in data_folders:
        base_name = f  # technically the name of the folder- doesn't need modification
        lengths = LengthTallyDict.copy()  # create copy of the template container for this specific file
    
        # Read data in, keep tallies
        bamfile = pysam.AlignmentFile("./Tophat_output/%s/%s/%s_sorted.bam" % (args.dataset, f, f), "rb")
        for read in bamfile.fetch(until_eof = True):
            lengths[read.query_length] += float(1)/read.get_tag("NH")
        bamfile.close()
    
    
        # Prepare first outstring for the counts file, write it out
        outstring = base_name
        for i in lengths.keys():
            outstring += "," + str(lengths[i])
        outfh_len.write(outstring + "\n")
    
    
        # Prepare second outstring for percents file, write it out
        lentot = sum(lengths.values())
        outstring = base_name
        for i in lengths.keys():
            outstring += "," + str(round(float(100*lengths[i])/lentot,3))
        outfh_lenperc.write(outstring + "\n")
    
        print "Finished %s at" % f, time.ctime()

    # Close output files when done with all samples
    outfh_len.close()
    outfh_lenperc.close()

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate the length distribution of the reads in a dataset after adapter trimming, rRNA+tRNA decontamination, and alignment.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Tophat_output folder (S1, S2, etc)")
    args = parser.parse_args()

    calculate_length_distribution(args.dataset)

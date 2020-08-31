#!/usr/bin/env python

from Bio import SeqIO
import argparse
import time
import collections
import sys
import os


def calculate_length_distribution(dataset):
    # Create a template container for the tallies of read lengths to be stored #
    LengthTallyDict = collections.OrderedDict()
    for i in range(0, 126):
        LengthTallyDict[i] = 0
    
    # Get list of input files, in fastq format #
    data_files = sorted(os.listdir("./Decontaminated_reads/%s/" % args.dataset))
    
    # Create folders for the output to be placed #
    if not os.path.exists("./output_length_dist_pre_selection"):
        os.mkdir("./output_length_dist_pre_selection")
    if not os.path.exists("./output_length_dist_pre_selection/%s" % args.dataset):
        os.mkdir("./output_length_dist_pre_selection/%s" % args.dataset)
        
    
    # Create output file (.csv) to place the counts of reads with each length
    # Header looks like Sample,18,19,20,21 ... 
    # First row looks like 1R_CAGATC_L002_R1_001,12345,23456, ... 
    outfh_len = open("./output_length_dist_pre_selection/%s/%s_pre_selection_length_counts.csv" % (args.dataset, args.dataset), "w")
    header_len = "Sample"
    for i in LengthTallyDict.keys():
        header_len += "," + str(i)
    outfh_len.write(header_len + "\n")


    # Create output file (.csv) to place the percents of reads with each length
    # Header looks like Sample,% Reads Len 18,% Reads Len 19, ... 
    # First row looks like 1R_CAGATC_L002_R1_001,12.345,23.456, ... 
    outfh_lenperc = open("./output_length_dist_pre_selection/%s/%s_pre_selection_length_percents.csv" % (args.dataset, args.dataset), "w")
    header_len = "Sample"
    for i in LengthTallyDict.keys():
        header_len += ",% Reads Len " + str(i)
    outfh_lenperc.write(header_len + "\n")
    
    
    for f in data_files:
        base_name = f.rsplit("_")[0]
        lengths = LengthTallyDict.copy()
    
        # Read data in, keep tallies
        infh = open("Decontaminated_reads/%s/%s" % (args.dataset, f), "r")
        for read in SeqIO.parse(infh, "fastq"):
            lengths[len(read)] += 1
        infh.close()
    
    
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
    
        print("Finished %s at" % f, time.ctime())

    # Close output files when done with all samples
    outfh_len.close()
    outfh_lenperc.close()

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate the length distribution of the reads in a dataset after adapter trimming/rRNA+tRNA decontamination, but before alignment and RPF length selection.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Decontaminated_reads folder (S1, S2, etc)")
    args = parser.parse_args()

    calculate_length_distribution(args.dataset)

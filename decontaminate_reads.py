#!/usr/bin/env python

import sys
import os
import math
import subprocess
import argparse
import time
from collections import OrderedDict


def create_output_folders(dataset):
    if not os.path.exists("./Decontaminated_reads"):
        os.mkdir("./Decontaminated_reads")
    if not os.path.exists("./Decontaminated_reads/%s" % dataset):
        os.mkdir("./Decontaminated_reads/%s" % dataset)
        
    if not os.path.exists("./Contaminated_read_alignments"):
        os.mkdir("./Contaminated_read_alignments")
    if not os.path.exists("./Contaminated_read_alignments/%s" % dataset):
        os.mkdir("./Contaminated_read_alignments/%s" % dataset)
        
    if not os.path.exists("./Logs_decontamination"):
        os.mkdir("./Logs_decontamination/")        
    if not os.path.exists("./Logs_decontamination/%s" % dataset):
        os.mkdir("./Logs_decontamination/%s" % dataset)


def unzip(dataset, base_name):
    unzip_cmd = "gunzip ./Trimmed/%s/%s_trimmed.fastq.gz" % (dataset, base_name)
    p2 = subprocess.Popen(unzip_cmd, shell = True)
    p2.wait()


def decontaminate_reads(dataset, base_name):
    print(dataset)
    print(base_name)
    bowtie_cmd = "bowtie -l 20 --un=Decontaminated_reads/%s/%s_decontaminated.fastq HisatIndex/contaminant_sequences Trimmed/%s/%s_trimmed.fastq >Contaminated_read_alignments/%s/%s_contaminated_read_alignments.aln" % (dataset, base_name, dataset, base_name, dataset, base_name)
    
    log_fh = open("./Logs_decontamination/%s/%s_decontamination.log" % (dataset, base_name), "w")
    print(bowtie_cmd)
    p = subprocess.Popen(bowtie_cmd, shell = True, stderr = log_fh)
    p.wait()

    log_fh.flush()
    log_fh.close()

    print("Finished decontaminating sample %s at" % base_name, time.ctime())
    print("Compressing %s..." % base_name)
    
    gzip_cmd = "gzip ./Trimmed/%s/%s_trimmed.fastq" % (dataset, base_name)
    p3 = subprocess.Popen(gzip_cmd, shell = True)
    p3.wait()
    
    print("Finished compressing %s at " % base_name, time.ctime())
    
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Decontaminate the rRNA/tRNA from adapter-trimmed reads.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Trimmed folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    
    data_files = sorted([x for x in os.listdir("./Trimmed/%s/" % args.dataset) if ".fastq" in x])
    for f in data_files:
        
        base_name = "_".join(f.split(".")[0].split("_")[:-1])
        if ".gz" in f:
            unzip(args.dataset, base_name)
            decontaminate_reads(args.dataset, base_name)
        else:
            decontaminate_reads(args.dataset, base_name)

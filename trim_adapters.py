#!/usr/bin/env python

import sys
import os
import math
import subprocess
import argparse
import time
from collections import OrderedDict


def create_output_folders(dataset):
    if not os.path.exists("./Trimmed"):
        os.mkdir("./Trimmed")
    if not os.path.exists("./Trimmed/%s" % dataset):
        os.mkdir("./Trimmed/%s" % dataset)
    if not os.path.exists("./Logs_trimmed"):
        os.mkdir("./Logs_trimmed")
    if not os.path.exists("./Logs_trimmed/%s" % dataset):
        os.mkdir("./Logs_trimmed/%s" % dataset)



def trim_adapters(dataset, base_name):
    
    zcat_cmd = "zcat Raw_data/%s/%s.fastq.gz" % (dataset, base_name)
    fastx_cmd = "fastx_clipper -a AGATCGGAAGAGCACACGTCT -l 18 -c -n -v -Q33 >Trimmed/%s/%s_trimmed.fastq" % (dataset, base_name)


    log_fh = open("./Logs_trimmed/%s/%s_trimmed.log" % (dataset, base_name), "w")

    p1 = subprocess.Popen(zcat_cmd, shell = True, stdout = subprocess.PIPE)
    p2 = subprocess.Popen(fastx_cmd, shell = True, stdin = p1.stdout, stderr = log_fh)
    p2.wait()

    log_fh.flush()
    log_fh.close()
    
    print("Finished trimming %s at " % base_name , time.ctime())
    print("Compressing %s..." % base_name)
    
    gzip_cmd = "gzip ./Trimmed/%s/%s_trimmed.fastq" % (dataset, base_name)
    p3 = subprocess.Popen(gzip_cmd, shell = True)
    p3.wait()
    
    print("Finished compressing %s at " % base_name, time.ctime())
    
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Trim the adapters from raw reads.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Raw_data folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    
    data_files = sorted([x for x in os.listdir("./Raw_data/%s/" % args.dataset) if ".fastq.gz" in x])
    if len(data_files) == 0:
        sys.exit("No files found in \'./Raw_data/%s/\'. Please ensure data files are present and in .fastq.gz format" % args.dataset)
    for f in data_files:
        base_name = f.split(".")[0]
        trim_adapters(args.dataset, base_name)

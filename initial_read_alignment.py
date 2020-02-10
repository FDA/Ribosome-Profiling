#!/usr/bin/env python

import sys
import os
import math
import subprocess
import argparse
import time
from collections import OrderedDict

class DataSetInfo(object):
    def __init__(self, genome, wttranscriptome, opttranscriptome1, opttranscriptome2, chx, wt, opt1, opt2, mismatch_allowance, max_multihits):
        self.genome = genome
        self.wttranscriptome = wttranscriptome
        self.opttranscriptome1 = opttranscriptome1
        self.opttranscriptome2 = opttranscriptome2
        self.chx = chx
        self.wt = wt  # list
        self.opt1 = opt1  # list
        self.opt2 = opt2  # list
        self.mismatch_allowance = mismatch_allowance
        self.num_alignments = max_multihits

def define_dataset_info():
    global datasetinfo
    datasetinfo = {}
    datasetinfo["S12"] = DataSetInfo("genome_plus_F9_constructs", "transcripts_mWTF9_pWTF9c", "transcripts_mWTF9_popt1F9c", None, False, ["1","2","3"], ["4","5","6"], [], "2", "20")
    return datasetinfo


def create_output_folders(dataset):
    if not os.path.exists("./Tophat_output"):
        os.mkdir("./Tophat_output")
    if not os.path.exists("./Tophat_output/%s" % dataset):
        os.mkdir("./Tophat_output/%s" % dataset)
        
        
        
def run_tophat(dataset, base_name):
    
   
    if base_name[0] in datasetinfo[dataset].wt:
        transcriptome = datasetinfo[dataset].wttranscriptome
    elif base_name[0] in datasetinfo[dataset].opt1:
        transcriptome = datasetinfo[dataset].opttranscriptome1
    elif base_name[0] in datasetinfo[dataset].opt2:
        transcriptome = datasetinfo[dataset].opttranscriptome2
    else:
        assert False  # shouldn't happen-- check to make sure DataSetInfo is populated correctly for this sample 
    
    #The following command may be changed depnding on the version of bowtie used to create the reference
    tophat_cmd = "tophat --bowtie1 -p 6 --no-novel-juncs -N %s -g %s -o ./Tophat_output/%s/%s ./BowtieIndex/%s ./Decontaminated_reads/%s/%s_decontaminated.fastq" % (datasetinfo[dataset].mismatch_allowance, datasetinfo[dataset].num_alignments, dataset, base_name, transcriptome, dataset, base_name)
    print tophat_cmd
    p1 = subprocess.Popen(tophat_cmd, shell=True)
    p1.wait()
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run tophat on adapter-trimmed, decontaminated reads.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Decontaminated_reads folder (S1, S2, etc)")
    #parser.add_argument("num_alignments", help = "The maximum number of alignments allowed for a read, option -g in tophat (default value = 20)", default = 20)
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    define_dataset_info()

    data_files = sorted([x for x in os.listdir("./Decontaminated_reads/%s/" % args.dataset) if ".fastq" in x])
    print data_files
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[:-1])
        run_tophat(args.dataset, base_name)

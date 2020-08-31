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
        self.wt = wt
        self.opt1 = opt1
        self.opt2 = opt2
        self.mismatch_allowance = mismatch_allowance
        self.num_alignments = max_multihits

def define_dataset_info():
    global datasetinfo
    datasetinfo = {}
    datasetinfo["S12"] = DataSetInfo("genome_plus_F9_constructs", "transcripts_mWTF9_pWTF9c", "transcripts_mWTF9_popt1F9c", None, False, ["1","2","3"], ["4","5","6"], [], "2", "20")
    return datasetinfo


def create_output_folders(dataset):
    if not os.path.exists("./Hisat_output"):
        os.mkdir("./Hisat_output")
    if not os.path.exists("./Hisat_output/%s" % dataset):
        os.mkdir("./Hisat_output/%s" % dataset)
        
        
        
def run_hisat(dataset, base_name):
    if not os.path.exists("./Hisat_output/%s/%s" % (dataset, base_name)):
        os.mkdir("./Hisat_output/%s/%s" % (dataset, base_name))

    if base_name[0] in datasetinfo[dataset].wt:
        transcriptome = datasetinfo[dataset].wttranscriptome
    elif base_name[0] in datasetinfo[dataset].opt1:
        transcriptome = datasetinfo[dataset].opttranscriptome1
    elif base_name[0] in datasetinfo[dataset].opt2:
        transcriptome = datasetinfo[dataset].opttranscriptome2
    else:
        # shouldn't happen-- check to make sure DataSetInfo is populated correctly for this sample
        assert False

    #Alignment command
    hisat_cmd = "hisat2 -p 6 -k %s -x ./HisatIndex/%s -U ./Decontaminated_reads/%s/%s_decontaminated.fastq | samtools view -Sbh > ./Hisat_output/%s/%s/accepted_hits.bam" % (datasetinfo[dataset].num_alignments, transcriptome, dataset, base_name, dataset, base_name)

    print(hisat_cmd)
    p1 = subprocess.Popen(hisat_cmd, shell=True)
    p1.wait()
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Run hisat2 on adapter-trimmed, decontaminated reads.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Decontaminated_reads folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    define_dataset_info()

    data_files = sorted([x for x in os.listdir("./Decontaminated_reads/%s/" % args.dataset) if ".fastq" in x])
    print(data_files)
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[:-1])
        run_hisat(args.dataset, base_name)

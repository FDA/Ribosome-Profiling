#!/usr/bin/env python

import sys
import os
import math
import subprocess
import argparse
import time
import pysam
from collections import OrderedDict


def samtools_sort_index(dataset, base_name):
    pysam.sort("-o", "./Tophat_output/%s/%s/%s_sorted.bam" % (dataset, base_name, base_name), "./Tophat_output/%s/%s/accepted_hits.bam" % (dataset, base_name))
    pysam.index("./Tophat_output/%s/%s/%s_sorted.bam" % (dataset, base_name, base_name))
    
    print "Finished with %s at " % base_name, time.ctime()
    
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Use samtools to sort and index the BAM files outputted by TopHat.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Tophat_output folder (S1, S2, etc)")
    args = parser.parse_args()

    data_folders = sorted(os.listdir("./Tophat_output/%s/" % args.dataset))  # fix
    for f in data_folders:
        base_name = f
        samtools_sort_index(args.dataset, base_name)
        
        
        

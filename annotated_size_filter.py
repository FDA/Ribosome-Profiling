#!/usr/bin/env python

import os
import argparse
import pysam


class DataSetSizeInfo(object):
    def __init__(self, sizes):
        # list of sizes
        self.sizes = sizes
        
        
def define_dataset_sizes():
    global allowed_sizes
    allowed_sizes = {}
    allowed_sizes["S12"] = {}
    allowed_sizes["S12"]["R"] = [20, 21, 22, 27, 28, 29, 30]
    allowed_sizes["S12"]["T"] = range(20, 126)
    
def create_output_folders(dataset):
    if not os.path.exists("./Annotated_size_filtered_reads"):
        os.mkdir("./Annotated_size_filtered_reads")
    if not os.path.exists("./Annotated_size_filtered_reads/%s" % dataset):
        os.mkdir("./Annotated_size_filtered_reads/%s" % dataset)
        
        
def size_filter_reads(dataset, base_name):
    sample_type = base_name[1]
    
    bamfile = pysam.AlignmentFile("./Annotated_reads/%s/annotated_%s.bam" % (dataset, base_name), "rb")
    outbam = pysam.AlignmentFile("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name), "wb", template = bamfile)
    
    for read in bamfile.fetch():
        if read.query_length in allowed_sizes[dataset][sample_type]:
            outbam.write(read)
            
    bamfile.close()
    outbam.close()
    
def index_new_bamfile(dataset, base_name):   
    pysam.index("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (args.dataset, base_name))
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Filter the annotated reads to the appropriate size for the dataset. Make sure the sizes are coded in define_dataset_sizes() before running!")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Annotated_reads folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    define_dataset_sizes()
    
    
    data_files = sorted(os.listdir("./Annotated_reads/%s/" % args.dataset))    
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[1:])
        size_filter_reads(args.dataset, base_name)
        index_new_bamfile(args.dataset, base_name)

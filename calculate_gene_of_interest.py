#!/usr/bin/env python

import pysam
import sys
import os
import time
import collections
import argparse as ap
import pdb
start_time = time.time()


def define_dataset_sizes():
    global allowed_sizes
    allowed_sizes["S12"] = {}
    allowed_sizes["S12"]["R"] = [20, 21, 22, 27, 28, 29, 30]
    allowed_sizes["S12"]["T"] = range(20, 126)

def ensembl_ID_converter(dataset):
    convert = {}
    
    with open("%s_ensembl_to_gene_id.tsv" % dataset, "r") as infh:
        for line in infh:
            s = line.rstrip().split()
            convert[s[0]] = s[1]
            
    return convert
    
def read_gene_lengths(dataset, convert):
    gene_lengths = {}
    
    with open("%s_ensembl_CDS_lengths.tsv" % dataset, "r") as infh:
        for line in infh:
            s = line.split()
            if "ENSTR" in s[0] or "ENSG" in s[0]:
                continue
            gene_lengths[convert[s[0]]] = int(s[1])
            
    return gene_lengths
    

def create_output_folders(dataset):
    if not os.path.exists("./output_goi"):
        os.mkdir("./output_goi")
    if not os.path.exists("./output_goi/%s" % dataset):
        os.mkdir("./output_goi/%s" % dataset)
        
        
def calculate_gene_of_interest_coverage(convert, gene_lengths, goi_id, dataset, base_name):
    sample_type = base_name[1]  # R or T
    read_totals = {}
    bamfile = pysam.AlignmentFile("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (args.dataset, base_name), "rb")
    read_totals = {}
    tot_infh = open("./Annotated_totals_per_read/%s/%s_read_counts.tsv" % (args.dataset, base_name), "r")
    tot_infh.next()  # skip header
    for line in tot_infh:
        s = line.rstrip().split("\t")
        read_totals[s[0]] = float(s[1])
    tot_infh.close()
    goi_counts = collections.OrderedDict()
    goi_count_categories = ["All"] + allowed_sizes[dataset][sample_type]
    for size in goi_count_categories:
        goi_counts[size] = collections.OrderedDict()
        for i in range(0, gene_lengths[goi_id]/3):
            goi_counts[size][i] = 0
             
    print goi_counts.keys()
        
    for read in bamfile.fetch(until_eof = True):
        length = read.query_length
        if length not in allowed_sizes[dataset][sample_type]:
            continue
        if read.reference_name[0:5] == "ENSTR" or read.reference_name[0:4] == "ENSG":
            continue
        if convert[str(read.reference_name)] != goi_id:
            continue
        if read.get_tag("al") == "UTR":
            continue

        read_nh = read_totals[read.qname]
        goi_counts["All"][read.get_tag("pp")] += (1/read_nh)
        goi_counts[length][read.get_tag("pp")] += (1/read_nh)
        
    bamfile.close()
    
    write_gene_of_interest(goi_id, goi_counts, dataset, base_name)
    
    
def write_gene_of_interest(goi_id, goi_counts, dataset, base_name):
    outfile = "./output_goi/%s/%s_%s_fractional_read_counts_by_codon.tsv" % (dataset, goi_id, base_name)
    print outfile
    outfh = open(outfile, "w")
    header = ["Codon number"]
    for key in goi_counts.keys():
        header.append(str(key))
    headerline = "\t".join(header) + "\n"
    outfh.write(headerline)
    for position in goi_counts["All"].keys():
        outline = str(position)
        for size in goi_counts.keys():
            outline += "\t" + str(round(goi_counts[size][position],2))
        outline += "\n"
        outfh.write(outline)
    
    outfh.close()   

if __name__ == "__main__":

    parser = ap.ArgumentParser(description = "Calculate the reads aligning to each codon of a gene of interest.")
    parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
    parser.add_argument("gene_of_interest", help = "The gene ID of the gene of interest. Can be either ensembl identifier (ex. ENST00000331789.5) or the gene identifier (ex. ACTB). [Note: the output file will always use the gene identifier]")
    args = parser.parse_args()
    
    define_dataset_sizes()
    convert = ensembl_ID_converter(args.dataset)
    

    if args.gene_of_interest in convert.keys():
        goi_id = convert[args.gene_of_interest]
    elif args.gene_of_interest in convert.values():
        goi_id = args.gene_of_interest
    else:
        print "Invalid gene of interest provided. Exiting program now."
        
    print goi_id
        

    gene_lengths = read_gene_lengths(args.dataset, convert)
    create_output_folders(args.dataset)
    data_files = sorted([x for x in os.listdir("./Annotated_size_filtered_reads/%s/" % args.dataset) if x[-4:] == ".bam"])
    

    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[0:-2])
        calculate_gene_of_interest_coverage(convert, gene_lengths, goi_id, args.dataset, base_name)

#!/usr/bin/env python

import pysam
from Bio import SeqIO
import sys
import os
import time
import collections
import argparse
start_time = time.time()

def ensembl_ID_converter(dataset):
    convert = {}
    
    with open("%s_ensembl_to_gene_id.tsv" % dataset, "r") as infh:
        for line in infh:
            s = line.rstrip().split()
            convert[s[0]] = s[1]
            
    return convert

def parse_sequence(dataset, ensembl_id):
    codons = []
    infh = open("./%s_all_CDS.fa" % dataset, "r")
    for record in SeqIO.parse(infh, "fasta"):
        if record.id == ensembl_id:
            # avoid partial CDS
            assert len(record.seq)%3 == 0
            for i in range(0,len(record.seq), 3):
                codons.append(str(record.seq[i:i+3]))
    
    # if false, codons is empty and you didn't find the sequence
    assert codons
    return codons
            


def parse_goi_files(dataset, data_files, base_names, goi_id):
    data = {}
    for base_name in base_names:
        data[base_name] = collections.OrderedDict()
        

    
    for base_name in base_names:
        infh = open("./output_goi/%s/%s_%s_fractional_read_counts_by_codon.tsv" % (dataset, goi_id, base_name), "r")
        header = infh.readline().rstrip()
        read_lengths = header.split("\t")[2:]

        for line in infh:
            s = line.rstrip().split("\t")
            # number, not letters
            codon_index = int(s[0])
            total = float(s[1])
            data[base_name][codon_index] = total
            
            
        infh.close()
                   
    return data
    
def calculate_averages(dataset, base_names, codons, data):
    averages = collections.OrderedDict()
    for base_name in base_names:
        averages[base_name] = float(sum(data[base_name].values()))/len(data[base_name].values())
        
    return averages
        
    
    
    
def write_formatted_goi_files(dataset, base_names, goi_id, codons, data, averages):
    outfh_ribo = open("./output_goi/%s/%s_%s_ribo_compiled.csv" % (dataset, dataset, goi_id), "w")
    outfh_total = open("./output_goi/%s/%s_%s_mRNA_compiled.csv" % (dataset, dataset, goi_id), "w")
    
    ribo_base_names = []
    total_base_names = []
    for base_name in base_names:
        if base_name[1] == "R":
            ribo_base_names.append(base_name)
        elif base_name[1] == "T":
            total_base_names.append(base_name)
        else:
            assert false
    
    ribo_outline = []
    for base_name in ribo_base_names:
        ribo_header = ["CODON","NUMBER",base_name.split("_")[0],"AVG","NORM", ""]
        ribo_outline += ribo_header
    total_outline = []
    for base_name in total_base_names:
        total_header = ["CODON","NUMBER",base_name.split("_")[0],"AVG","NORM", ""]
        total_outline += total_header
        
    outfh_ribo.write(",".join(ribo_outline) + "\n")
    outfh_total.write(",".join(total_outline) + "\n")
    
    for i in range(0, len(codons)):
        ribo_outline = []
        for base_name in ribo_base_names:
            try:
                ribo_block = [str(codons[i]), str(i), str(data[base_name][i]), str(averages[base_name]), str(round(float(data[base_name][i])/averages[base_name],3)), ""]
                ribo_outline += ribo_block
            except ZeroDivisionError as e:
                ribo_block = [str(codons[i]), str(i), str(data[base_name][i]), str(averages[base_name]), str(0), ""]
                ribo_outline += ribo_block
            
        total_outline = []
        for base_name in total_base_names:
            try:
                total_block = [str(codons[i]), str(i), str(data[base_name][i]), str(averages[base_name]), str(round(float(data[base_name][i])/averages[base_name],3)), ""]
                total_outline += total_block
            except ZeroDivisionError as e:
                total_block = [str(codons[i]), str(i), str(data[base_name][i]), str(averages[base_name]), str(0), ""]
                total_outline += total_block
            
        outfh_ribo.write(",".join(ribo_outline) + "\n")
        outfh_total.write(",".join(total_outline) + "\n")
        
    outfh_ribo.close()
    outfh_total.close()
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Format gene of interest files into Gaya's format and do some calculations.")
    parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
    parser.add_argument("gene_of_interest", help = "The gene of interest. There MUST already be files in ./output_goi/S#/ for this gene")
    args = parser.parse_args()
    
    convert = ensembl_ID_converter(args.dataset)
    
    if args.gene_of_interest in convert.keys():
        ensembl_id = args.gene_of_interest
        goi_id = convert[args.gene_of_interest]
    elif args.gene_of_interest in convert.values():
        goi_id = args.gene_of_interest
        for key in convert.keys():
            if convert[key] == goi_id:
                ensembl_id = key
    else:
        print("Invalid gene of interest provided. Exiting program now.")
        
    print(goi_id)
    print(ensembl_id)

    codons = parse_sequence(args.dataset, ensembl_id)
    
    data_files = sorted([x for x in os.listdir("./output_goi/%s/" % args.dataset) if x.split("_")[0] == args.gene_of_interest])
    
    base_names = []
    for f in data_files:
        base_name = f.split(".")[0].split("_")[1:-5]
        base_names.append("_".join(base_name))
        
        
    data = parse_goi_files(args.dataset, data_files, base_names, goi_id)
    averages = calculate_averages(args.dataset, base_names, codons, data)
    write_formatted_goi_files(args.dataset, base_names, goi_id, codons, data, averages)

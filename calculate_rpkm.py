#!/usr/bin/env python

import sys
import os
import collections
import argparse as ap
import pysam
import time
import math
import pdb
start_time = time.time()

ProteinTallyDict = collections.OrderedDict()

class SampleRPKM(object):
    def __init__(self, cutoff):
        self.protein_counts = ProteinTallyDict.copy()
        self.rpkm = ProteinTallyDict.copy()
        self.total_reads = 0
        self.cutoff = cutoff
        
    def write_output(self, outfh, cutoff, convert):
        outfh.write(",".join(["Gene", "GeneID", "RPKM", "Length (codons) (minus ends)", "# Reads in Gene", "Total Reads in Sample"]) + "\n")
        
        for protein in self.protein_counts.keys():
            try:
                if protein_lengths[protein]-(2*cutoff) <= 0:
                    outfh.write(",".join([protein, convert[protein], str(0), str(0), str(0), str(self.total_reads)]))
                else:
                    outfh.write(",".join([protein, convert[protein], str(self.rpkm[protein]), str(protein_lengths[protein]-(2*cutoff)), str(self.protein_counts[protein]), str(self.total_reads)]) + "\n")
            except KeyError as e:
                continue # because of transcripts without cdses


def create_output_folders(dataset):
    if not os.path.exists("./output_rpkm"):
        os.mkdir("./output_rpkm")
    if not os.path.exists("./output_rpkm/%s" % dataset):
        os.mkdir("./output_rpkm/%s" % dataset)
        
        
def get_start_end_coords(dataset):
    coordinates_dict = {}  # key: transcript id, value: dict w/ keys "start" and "end"
    infh = open("%s_transcript_CDS_start_end_locations.txt" % dataset, "r")
    infh.next()  # skip header
    for line in infh:
        s = line.split("\t")
        coordinates_dict[s[0]] = {"start":int(s[1])-1, "end":int(s[2])}  # subtract 1 from start bc the values are 1-indexed, end is inclusive of last base so it balances the -1 shift
                
    return coordinates_dict
    
    
def ensembl_ID_converter(dataset):
    convert = {}
    
    with open("%s_ensembl_to_gene_id.tsv" % dataset, "r") as infh:
        for line in infh:
            s = line.rstrip().split("\t")
            convert[s[0]] = s[1]  # key= ensembl_id, value= gene_id
            
    return convert
    
    
def read_transcript_lengths(convert, dataset):
    transcript_lengths = {}
    
    with open("%s_ensembl_transcript_lengths.tsv" % dataset, "r") as infh:
        for line in infh:
            s = line.split()

            transcript_lengths[s[0]] = int(s[1])

    return transcript_lengths 
    
    
def get_protein_lengths(coordinates_dict):
    protein_lengths = {}

    for key in coordinates_dict.keys():
        protein_lengths[key] = int(math.ceil((coordinates_dict[key]["end"] - coordinates_dict[key]["start"])/3))
    return protein_lengths 
    
    
def calculate_rpkm(dataset, base_name, sample_data, cutoffs, protein_lengths):
    sample_data[base_name] = {}
    for cutoff in cutoffs:
        sample_data[base_name][cutoff] = SampleRPKM(cutoff)

    read_totals = {}
    tot_infh = open("./Annotated_totals_per_read/%s/%s_read_counts.tsv" % (dataset, base_name), "r")
    tot_infh.next()  # skip header
    for line in tot_infh:
        s = line.rstrip().split("\t")
        read_totals[s[0]] = float(s[1])
    tot_infh.close()    
    
    bamfile = pysam.AlignmentFile("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name), "rb")
    
    for read in bamfile.fetch(until_eof = True):
        A_site_pos = int(read.get_tag("pp"))
        
        if read.reference_name[0:5] == "ENSTR" or read.reference_name[0:4] == "ENSG":
            continue
        
        for cutoff in cutoffs:
            if (A_site_pos < cutoff) or (A_site_pos > protein_lengths[read.reference_name]-cutoff):
                pass
            else:
                sample_data[base_name][cutoff].protein_counts[read.reference_name] += 1/read_totals[read.qname]

    for cutoff in cutoffs:
        sample_data[base_name][cutoff].total_reads = sum(sample_data[base_name][cutoff].protein_counts.values())
    
        for protein in sample_data[base_name][cutoff].protein_counts.keys():
        
        
            if protein[0:5] == "ENSTR" or protein[0:4] == "ENSG":
                continue
                
            try:    
            
            
                if (protein_lengths[protein]-(cutoff*2)) <= 0:
                    sample_data[base_name][cutoff].rpkm[protein] = 0
                else:
                    sample_data[base_name][cutoff].rpkm[protein] = round(float(1000000000 * sample_data[base_name][cutoff].protein_counts[protein])/(sample_data[base_name][cutoff].total_reads * (protein_lengths[protein]-(cutoff*2))),4)
                    
                
            except KeyError as e:
                continue   
                # This exists because there are transcripts without CDSs in the alignment
                
    return sample_data

if __name__ == "__main__":
    parser = ap.ArgumentParser(description = "Calculate the RPKM.")
    parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    convert = ensembl_ID_converter(args.dataset)
    transcript_lengths = read_transcript_lengths(convert, args.dataset)
    coordinates_dict = get_start_end_coords(args.dataset)
    protein_lengths = get_protein_lengths(coordinates_dict)
    
    for key in transcript_lengths.keys():
        ProteinTallyDict[key] = 0
        
    cutoffs = [0, 20, 30, 40, 50]
    sample_data = {}
    base_names = []
    
    data_files = sorted(os.listdir("./Annotated_size_filtered_reads/%s/" % args.dataset))
    
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[0:-2])

        base_names.append(base_name)
        sample_data = calculate_rpkm(args.dataset, base_name, sample_data, cutoffs, protein_lengths)
        
        
        
        
    for base_name in base_names:
        for cutoff in cutoffs:
            outfh = open("./output_rpkm/%s/%s_RPKM_cutoff_%s.csv" % (args.dataset, base_name, str(cutoff)), "w")
            sample_data[base_name][cutoff].write_output(outfh, cutoff, convert)
            outfh.close()
        
        
    for cutoff in cutoffs:
        file_columns = ["RPKM %s" % base_name for base_name in base_names]
        header = "Gene,GeneID," + ",".join(file_columns) + "\n"
        outfh = open("./output_rpkm/%s/%s_RPKM_comparison_cutoff_%s.csv" % (args.dataset, args.dataset, str(cutoff)), "w")
        outfh.write(header)
        
        for gene in ProteinTallyDict.keys():
            output_list = [gene, convert[gene]]
            for base_name in base_names:
                output_list.append(str(sample_data[base_name][cutoff].rpkm[gene]))
            
            outfh.write(",".join(output_list) + "\n")
            
        outfh.close()


#/usr/bin/env python

import pysam
import sys
import os
import time
import collections
from Bio import SeqIO
import argparse as ap
import pdb
from initial_read_alignment import DataSetInfo
from initial_read_alignment import define_dataset_info
from annotated_size_filter import define_dataset_sizes
start_time = time.time()


special_cases = {}
special_cases["S12"] = ["F9WT", "F9opt1"]

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
    
def read_transcript_lengths(dataset, convert):
    transcript_lengths = {}
    
    for record in SeqIO.parse(open("%s_all_transcripts.fa" % dataset, "r"), "fasta"):
        #print record.id, record.description
        if "ENSTR" in record.id or "ENSG" in record.id:
            continue
        transcript_lengths[convert[record.id]] = len(record.seq)
        
    return transcript_lengths
    
    
def get_start_end_coords(dataset):
    # coordinates_dict is... key: transcript id, value: dict w/ keys "start" and "end"
    coordinates_dict = {}
    infh = open("%s_transcript_CDS_start_end_locations.txt" % dataset, "r")
    next(infh)
    for line in infh:
        s = line.split("\t")
        # subtract 1 from start bc the values are 1-indexed, end is inclusive of last base so it balances the -1 shift
        coordinates_dict[s[0]] = {"start":int(s[1])-1, "end":int(s[2])}
                
    return coordinates_dict
    
    
def get_mRNA_sequence(dataset, goi_id, transcript_lengths):
    length = transcript_lengths[convert[goi_id]]
    sequence_dict = {}
    for i in range(0, length):
        sequence_dict[i] = "X"

    for record in SeqIO.parse(open("%s_all_transcripts.fa" % dataset, "r"), "fasta"):
        if record.id == goi_id:
            for i in range(0, len(record.seq)):
                sequence_dict[i] = record.seq[i]
                
                
    return sequence_dict
    
def get_CDS_sequence(dataset, goi_id, gene_lengths):
    length = gene_lengths[convert[goi_id]]
    cds_sequence_dict = {}
    for i in range(0, length):
        cds_sequence_dict[i] = "X"
        
    for record in SeqIO.parse(open("%s_all_CDS.fa" % dataset, "r"), "fasta"):
        if record.id == goi_id:
            for i in range(0, len(record.seq)):
                cds_sequence_dict[i] = record.seq[i]
                
    return cds_sequence_dict
    

def create_output_folders(dataset):
    if not os.path.exists("./output_mRNA_coverage"):
        os.mkdir("./output_mRNA_coverage")
    if not os.path.exists("./output_mRNA_coverage/%s" % dataset):
        os.mkdir("./output_mRNA_coverage/%s" % dataset)
        
        
def calculate_mRNA_coverage(convert, transcript_lengths, coordinates_dict, sequence_dict, goi_id, dataset, base_name):

    bamfile = pysam.AlignmentFile("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (args.dataset, base_name), "rb")

    
    if convert[goi_id] in special_cases[dataset]:
        print("Special")
        print(goi_id)
        print(convert[goi_id])
        if (convert[goi_id] == "ADAMTS13WT") or (convert[goi_id] == "F9WT"):
            print(base_name, convert[goi_id])
            
            if base_name[0] in datasetinfo[dataset].wt:
                print("yep it's here")
                string = goi_id + ":" + str(0) + "-" +  str(transcript_lengths[convert[goi_id]])
                mpileup = pysam.mpileup("-a", "-r", string, "./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name)).split("\n")            
                mrna_coverage = {}
                
                if len(mpileup) == 1 and mpileup[0] == "":
                    for i in range(0, transcript_lengths[convert[goi_id]]):
                        mrna_coverage[i] = 0
                    print("Finished filling in coverage with zeroes because no reads aligned to this gene for %s"  % base_name)
                    
                for i in range(0, len(mpileup)):
                    if mpileup[i] == "":
                        continue
                    s = mpileup[i].split("\t")
                    try:
                        mrna_position = int(s[1]) -1
                    except:
                        print(s)
                        print(base_name)
                    coverage = int(s[3])
                    mrna_coverage[i] = int(s[3])
                    
                print("Finished calculating mRNA coverage for %s" % base_name, "at", time.ctime())
                return mrna_coverage
            else:
                mrna_coverage = {}
                for i in range(0, transcript_lengths[convert[goi_id]]):
                    mrna_coverage[i] = 0
                print("Finished fudging mRNA coverage for %s" % base_name, "at", time.ctime())
                return mrna_coverage
                
        elif (convert[goi_id] == "ADAMTS13P118P") or (convert[goi_id] == "F9opt1") or (convert[goi_id] == "F91A") or (convert[goi_id] == "F9V107V"): 
            print(base_name, convert[goi_id])
            if base_name[0] in datasetinfo[dataset].opt1:
                print("yep it's here")
                string = goi_id + ":" + str(0) + "-" +  str(transcript_lengths[convert[goi_id]])
                mpileup = pysam.mpileup("-a", "-r", string, "./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name)).split("\n")            
                mrna_coverage = {}
                
                if len(mpileup) == 1 and mpileup[0] == "":
                    for i in range(0, transcript_lengths[convert[goi_id]]):
                        mrna_coverage[i] = 0
                    print("Finished filling in coverage with zeroes because no reads aligned to this gene for %s"  % base_name)
                    
                for i in range(0, len(mpileup)):
                    if mpileup[i] == "":
                        continue
                    s = mpileup[i].split("\t")
                    try:
                        mrna_position = int(s[1]) -1
                    except:
                        print(s)
                        print(base_name)
                    coverage = int(s[3])
                    mrna_coverage[i] = int(s[3])
                    
                print("Finished calculating mRNA coverage for %s" % base_name, "at", time.ctime())
                return mrna_coverage
            else:
                mrna_coverage = {}
                for i in range(0, transcript_lengths[convert[goi_id]]):
                    mrna_coverage[i] = 0
                print("Finished fudging mRNA coverage for %s" % base_name, "at", time.ctime())
                return mrna_coverage
                
        elif (convert[goi_id] == "F92A"):
            print("hello")
            print(base_name, convert[goi_id])
            print(datasetinfo[dataset].opt2)
            if base_name[0] in datasetinfo[dataset].opt2:
                string = goi_id + ":" + str(0) + "-" +  str(transcript_lengths[convert[goi_id]])
                mpileup = pysam.mpileup("-a", "-r", string, "./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name)).split("\n") 
                mrna_coverage = {}

                if len(mpileup) == 1 and mpileup[0] == "":
                    for i in range(0, transcript_lengths[convert[goi_id]]):
                        mrna_coverage[i] = 0
                    print("Finished filling in coverage with zeroes because no reads aligned to this gene for %s"  % base_name )

                for i in range(0, len(mpileup)):
                    if mpileup[i] == "":
                        continue
                    s = mpileup[i].split("\t")
                    try:
                        mrna_position = int(s[1]) -1
                    except:
                        print(s)
                        print(base_name)
                    coverage = int(s[3])
                    mrna_coverage[i] = int(s[3])
                    
                print("Finished calculating mRNA coverage for %s" % base_name, "at", time.ctime())
                return mrna_coverage
            
            else:
                mrna_coverage = {}
                for i in range(0, transcript_lengths[convert[goi_id]]):
                    mrna_coverage[i] = 0
                print("Finished fudging mRNA coverage for %s" % base_name, "at", time.ctime())
                return mrna_coverage                

            
    else:
    
        string = goi_id + ":" + str(0) + "-" +  str(transcript_lengths[convert[goi_id]])
        mpileup = pysam.mpileup("-a", "-r", string, "./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name)).split("\n")
        
        mrna_coverage = {}
        
        print(base_name)

        if len(mpileup) == 1 and mpileup[0] == "":
            for i in range(0, transcript_lengths[convert[goi_id]]):
                mrna_coverage[i] = 0
            print("Finished filling in coverage with zeroes because no reads aligned to this gene for %s"  % base_name)

        else:
        
            for i in range(0, len(mpileup)):
                if mpileup[i] == "":
                    continue
                s = mpileup[i].split("\t")
                try:
                    mrna_position = int(s[1]) -1
                except:
                    print(s)
                    print(base_name)
                coverage = int(s[3])
                mrna_coverage[i] = int(s[3])
                
                
            print("Finished calculating mRNA coverage for %s" % base_name, "at", time.ctime())

        return mrna_coverage
       
       
def write_mrna_coverage(mrna_coverage_dict, sequence_dict, coordinates_dict, transcript_lengths, convert, goi_id, dataset, base_names):
    
    outfh = open("./output_mRNA_coverage/%s/%s_%s_mRNA_coverage.csv" % (dataset, dataset, convert[goi_id]), "w")
    headers = ["mRNA Position", "CDS position", "Nucleotide"]
    for base_name in base_names:
        headers.append(base_name)
    outfh.write(",".join(headers) + "\n")
    
    for i in range(0, transcript_lengths[convert[goi_id]]):
        mrna_position = i
        cds_position = mrna_position - coordinates_dict[goi_id]["start"]
        nucleotide = sequence_dict[i]
        outline = [mrna_position, cds_position, nucleotide]
        for base_name in base_names:
            outline.append(mrna_coverage_dict[base_name][i])
            
        outstring = ",".join([str(x) for x in outline]) + "\n"
        outfh.write(outstring)
        
    outfh.close()
            
        
        
        
def write_mrna_coverage_codon_averages(mrna_coverage_dict, cds_sequence_dict, coordinates_dict, gene_lengths, convert, goi_id, dataset, base_names):

    outfh = open("./output_mRNA_coverage/%s/%s_%s_mRNA_coverage_codon_averages.csv" % (dataset, dataset, convert[goi_id]), "w")
    headers = ["Codon Number", "Codon"]
    for base_name in base_names:
        headers.append(base_name)
    outfh.write(",".join(headers) + "\n")
    
    for i in range(0, gene_lengths[convert[goi_id]], 3):
        codon_number = i/3
        codon = cds_sequence_dict[i] + cds_sequence_dict[i+1] + cds_sequence_dict[i+2]
        outline = [codon_number, codon]
        for base_name in base_names:
            coverage_values = [mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]],mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]+1],mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]+2]]
            coverage_average = round(float(sum(coverage_values))/3,2)
            outline.append(coverage_average)
            
        outstring = ",".join([str(x) for x in outline]) + "\n"
        outfh.write(outstring)
        
    outfh.close()
    
    
    
def write_normalized_mrna_coverage_codon_averages(mrna_coverage_dict, cds_sequence_dict, coordinates_dict, gene_lengths, convert, goi_id, dataset, base_names):

    all_codon_averages = {}
    for base_name in base_names:
        all_codon_averages[base_name] = []

    outfh = open("./output_mRNA_coverage/%s/%s_%s_normalized_mRNA_coverage_codon_averages.csv" % (dataset, dataset, convert[goi_id]), "w")
    headers = ["Codon Number", "Codon"]
    for base_name in base_names:
        headers.append(base_name)
    outfh.write(",".join(headers) + "\n")
    
    for i in range(0, gene_lengths[convert[goi_id]], 3):
        codon_number = i/3
        codon = cds_sequence_dict[i] + cds_sequence_dict[i+1] + cds_sequence_dict[i+2]
        outline = [codon_number, codon]
        for base_name in base_names:
            coverage_values = [mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]],mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]+1],mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]+2]]
            coverage_average = float(sum(coverage_values))/3
            all_codon_averages[base_name].append(coverage_average)
            

    normalization_numbers = {}
    for base_name in base_names:
        normalization_numbers[base_name] = float(sum(all_codon_averages[base_name]))/len(all_codon_averages[base_name])
        print(base_name, normalization_numbers[base_name])
    
    for i in range(0, gene_lengths[convert[goi_id]], 3):
        codon_number = i/3
        codon = cds_sequence_dict[i] + cds_sequence_dict[i+1] + cds_sequence_dict[i+2]
        outline = [codon_number, codon]
        for base_name in base_names:
            coverage_values = [mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]],mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]+1],mrna_coverage_dict[base_name][i+coordinates_dict[goi_id]["start"]+2]]
            coverage_average = float(sum(coverage_values))/3
            try:
                out_value = round(float(coverage_average)/normalization_numbers[base_name],5)
            except ZeroDivisionError as e:
                out_value = 0
            outline.append(str(out_value))
            
        outstring = ",".join([str(x) for x in outline]) + "\n"
        outfh.write(outstring)
        
    outfh.close()
        
  

if __name__ == "__main__":

    parser = ap.ArgumentParser(description = "Calculate the reads aligning to each codon of a gene of interest.")
    parser.add_argument("dataset", help = "The dataset to be looked at, e.g. S1")
    parser.add_argument("gene_of_interest", help = "The gene ID of the gene of interest. Can be either ensembl identifier (ex. ENST00000331789.5) or the gene identifier (ex. ACTB). [Note: the output file will always use the gene identifier]")
    args = parser.parse_args()
    
    global datasetinfo
    datasetinfo = define_dataset_info()
    define_dataset_sizes()
    convert = ensembl_ID_converter(args.dataset)
    

    if args.gene_of_interest in convert.keys():
        goi_id = args.gene_of_interest
    elif args.gene_of_interest in convert.values():
        for key in convert.keys():
            if convert[key] == args.gene_of_interest:
                goi_id = key
    else:
        print("Invalid gene of interest provided. Exiting program now.")

        
    print(goi_id)
        

    transcript_lengths = read_transcript_lengths(args.dataset, convert)
    gene_lengths = read_gene_lengths(args.dataset, convert)
    coordinates_dict = get_start_end_coords(args.dataset)
    sequence_dict = get_mRNA_sequence(args.dataset, goi_id, transcript_lengths)
    cds_sequence_dict = get_CDS_sequence(args.dataset, goi_id, gene_lengths)
    create_output_folders(args.dataset)
    
    data_files = sorted([x for x in os.listdir("./Annotated_size_filtered_reads/%s" % args.dataset) if x[-4:] == ".bam"])  # fix
    
    print(data_files)
    
    mrna_coverage_dict = {}
    base_names = []
    

    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[0:-2])
        base_names.append(base_name)
        mrna_coverage_dict[base_name] = calculate_mRNA_coverage(convert, transcript_lengths, coordinates_dict, sequence_dict, goi_id, args.dataset, base_name)        
    
    write_mrna_coverage(mrna_coverage_dict, sequence_dict, coordinates_dict, transcript_lengths, convert, goi_id, args.dataset, base_names)
    write_mrna_coverage_codon_averages(mrna_coverage_dict, cds_sequence_dict, coordinates_dict, gene_lengths, convert, goi_id, args.dataset, base_names)
    write_normalized_mrna_coverage_codon_averages(mrna_coverage_dict, cds_sequence_dict, coordinates_dict, gene_lengths, convert, goi_id, args.dataset, base_names)

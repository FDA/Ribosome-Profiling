#!/usr/bin/env python

import sys
import os
import collections
import argparse
from Bio import SeqIO
import pysam
import time
import math
import pdb
start_time = time.time()

                          
CodonsDict = collections.OrderedDict([('ATG', 0), ('TTT', 0), ('TTC', 0), ('TTA', 0), ('TTG', 0),
                          ('CTT', 0), ('CTC', 0), ('CTA', 0), ('CTG', 0), ('ATT', 0), ('ATC', 0),
                          ('ATA', 0), ('GTT', 0), ('GTC', 0), ('GTA', 0), ('GTG', 0), ('TCT', 0),
                          ('TCC', 0), ('TCA', 0), ('TCG', 0), ('AGT', 0), ('AGC', 0), ('CCT', 0),
                          ('CCC', 0), ('CCA', 0), ('CCG', 0), ('ACT', 0), ('ACC', 0), ('ACA', 0),
                          ('ACG', 0), ('GCT', 0), ('GCC', 0), ('GCA', 0), ('GCG', 0), ('TAT', 0),
                          ('TAC', 0), ('CAT', 0), ('CAC', 0), ('AAT', 0), ('AAC', 0), ('AAA', 0),
                          ('AAG', 0), ('GAT', 0), ('GAC', 0), ('GAA', 0), ('GAG', 0), ('TGT', 0),
                          ('TGC', 0), ('TGG', 0), ('CGT', 0), ('CGC', 0), ('CGA', 0), ('CGG', 0),
                          ('AGA', 0), ('AGG', 0), ('GGT', 0), ('GGC', 0), ('GGA', 0), ('GGG', 0),
                          ('CAA', 0), ('CAG', 0), ('TAG', 0), ('TGA', 0), ('TAA', 0)]) # initialize w/ acceptable keys
                          
                          
def create_output_folders(dataset):
    if not os.path.exists("./output_global_codon_patterns"):
        os.mkdir("./output_global_codon_patterns")
    if not os.path.exists("./output_global_codon_patterns/%s" % dataset):
        os.mkdir("./output_global_codon_patterns/%s" % dataset)
        
    if not os.path.exists("./output_global_codon_patterns_nondivided"):
        os.mkdir("./output_global_codon_patterns_nondivided")
    if not os.path.exists("./output_global_codon_patterns_nondivided/%s" % dataset):
        os.mkdir("./output_global_codon_patterns_nondivided/%s" % dataset)
        
        
        
        
def read_proteins(protein_filepath):  # this is needed because didn't use only complete CDSs
    protein_sequences = {}  # key = transcript_id, value = protein sequence
    protein_lengths = {}  # key = transcript_id, value = protein length (which == CDS length in codons)
    
    infh = open(protein_filepath, "r")
    for record in SeqIO.parse(infh, "fasta"):
        protein_sequences[record.id] = record.seq
        protein_lengths[record.id] = len(record.seq)
        
    return protein_sequences, protein_lengths


def calculate_codon_usage(cds_filepath, protein_sequences):
    codon_usage = {}  # key = transcript id, value = CodonsDict copy with codons filled in
    
    infh = open(cds_filepath, "r")
    for record in SeqIO.parse(infh, "fasta"):
        transcript_id = record.id  
        codon_usage[transcript_id] = CodonsDict.copy()
        
        if len(record.seq)%3 != 0:  # check to make sure CDS is complete
            for i in range(0,3):
                if record.seq[i:].translate(stop_symbol=".") == protein_sequences[record.id]:
                    frame = i
                    #SeqIO.write(record, outfh, "fasta")
        else:
            frame = 0 
            try:
                assert record.seq[0:].translate(stop_symbol=".") == protein_sequences[record.id]
            except AssertionError as e:
                print record.seq[0:].translate()
                print record.seq[0:].translate(stop_symbol=".")
                print protein_sequences[record.id]
                raise
                
        
        for i in range(frame, len(record.seq), 3):
            codon = record.seq[i:i+3].upper()
            
            if len(codon)<3:
                break
            
            codon_usage[transcript_id][codon] += 1

    return codon_usage
    
    
def get_start_end_coords(dataset):
    coordinates_dict = {}  # key: transcript id, value: dict w/ keys "start" and "end"
    infh = open("%s_transcript_CDS_start_end_locations.txt" % dataset, "r")
    infh.next()  # skip header
    for line in infh:
        s = line.split("\t")
        coordinates_dict[s[0]] = {"start":int(s[1])-1, "end":int(s[2])}  # subtract 1 from start bc the values are 1-indexed, end is inclusive of last base so it balances the -1 shift
        
    return coordinates_dict
    
    
def ensembl_ID_converter(dataset):
    convert = {}  # key= ensembl transcript_id, value= gene_id
    
    with open("%s_ensembl_to_gene_id.tsv" % dataset, "r") as infh:
        for line in infh:
            s = line.rstrip().split()
            convert[s[0]] = s[1]  
            
    return convert
    
def read_rpkm_values(dataset, base_name, cutoff):
    sample_rpkms = {}  # key = transcript_id, value = rpkm (float)
    
    infh = open("./output_rpkm/%s/%s_RPKM_cutoff_%s.csv" % (dataset, base_name, str(cutoff)), "r")
    infh.next()  # skip header line
    
    for line in infh:
        s = line.split(",")
        sample_rpkms[s[0]] = float(s[2])

    return sample_rpkms

def calculate_global_codon_patterns(dataset, base_name, sample_rpkms, codon_usage, protein_lengths, convert, averagelines):
    sites = ["A", "P", "E", "F"]
    sample_read_codons = {}
    sample_read_totals = {}  # total reads per gene    
    for site in sites:
        sample_read_codons[site] = {}
        
    for site in sites:
        for key in sample_rpkms.keys():
            sample_read_codons[site][key] = CodonsDict.copy()
            sample_read_totals[key] = 0
        
    read_totals = {}
    tot_infh = open("./Annotated_totals_per_read/%s/%s_read_counts.tsv" % (dataset, base_name), "r")
    tot_infh.next()  # skip header
    for line in tot_infh:
        s = line.rstrip().split("\t")
        read_totals[s[0]] = float(s[1])
    tot_infh.close()  
        
        
    bamfile = pysam.AlignmentFile("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name), "rb")
    for read in bamfile.fetch(until_eof = True):
        if read.reference_name[0:5] == "ENSTR" or read.reference_name[0:4] == "ENSG":
            continue
        read_gn = read.get_tag("gn")
        read_ac = read.get_tag("ac")  # A-site codon
        read_pc = read.get_tag("pc")  # P-site codon
        read_ec = read.get_tag("ec")  # E-site codon
        try:
            read_fc = read.get_tag("fc")  # "far"-site codon
        except:
            continue
        read_al = read.get_tag("al")  # al for "alignment type", CDS or UTR
        
        if read_al == "CDS":
            sample_read_codons["A"][read_gn][read_ac] += float(1)/read_totals[read.query_name]
            sample_read_codons["P"][read_gn][read_pc] += float(1)/read_totals[read.query_name]
            sample_read_codons["E"][read_gn][read_ec] += float(1)/read_totals[read.query_name]
            try:
                sample_read_codons["F"][read_gn][read_fc] += float(1)/read_totals[read.query_name]
            except:
                continue
            sample_read_totals[read_gn] += float(1)/read_totals[read.query_name]
            
    codon_values = {}
    for site in sites:
        codon_values[site] = {}
    for site in sites:
        for key in CodonsDict.keys():
            codon_values[site][key] = []
            
    for site in sites:
        total_reads_in_high_rpkm_cds = 0
        total_length_high_rpkm_cds = 0
        for key in sample_rpkms.keys():
            if sample_rpkms[key] < 10:
                continue
                
            total_reads_to_CDS = sample_read_totals[key]
            length_of_CDS = protein_lengths[key]    
            total_reads_in_high_rpkm_cds += sample_read_totals[key]
            total_length_high_rpkm_cds += protein_lengths[key]
            
            for c in sample_read_codons[site][key].keys():
                codon_c_in_gene = codon_usage[key][c]
                site_with_codon_c = sample_read_codons[site][key][c]
                
                try:
                    normalized_codon_freq = (site_with_codon_c)/((total_reads_to_CDS/length_of_CDS) * codon_c_in_gene)
                    codon_values[site][c].append(normalized_codon_freq)
                except ZeroDivisionError as e:
                    if (site_with_codon_c == total_reads_to_CDS) and (site_with_codon_c == 0):
                        codon_values[site][c].append(0)
                    elif (site_with_codon_c == 0) and (codon_c_in_gene == 0):
                        #outlist.append("K")
                        codon_values[site][c].append(0)
                    elif (site_with_codon_c > 0) and (codon_c_in_gene == 0):
                        #outlist.append("M")
                        codon_values[site][c].append(0)
                    elif (site_with_codon_c > 0) and (total_reads_to_CDS == 0):
                        #outlist.append("B")
                        codon_values[site][c].append(0)
                    
                    else:
                        codon_values[site][c].append("X")
        
        
    for site in sites:    
        outfh = open("./output_global_codon_patterns/%s/%s_%s_site_codon_patterns.csv" % (dataset, base_name, site), "w")
        headerlist = ["Gene", "Gene ID", "RPKM", "Total Reads in CDS", "Length of CDS"]
        for codon in CodonsDict.keys():
            headerlist.append("%s" % codon.replace("T","U"))
        outfh.write(",".join(headerlist) + "\n")
        
        averageline = ["All", "All", ">=10", str(total_reads_in_high_rpkm_cds), str(total_length_high_rpkm_cds)]
        for codon in CodonsDict.keys():
            averageline.append(str(sum(codon_values[site][codon])/len(codon_values[site][codon])))
        outfh.write(",".join(averageline) + "\n")
        averagelines[site][base_name] = averageline
        
        for key in sample_rpkms.keys():
            total_reads_to_CDS = sample_read_totals[key]
            length_of_CDS = protein_lengths[key]
            
            outlist = [str(key), str(convert[key]), str(sample_rpkms[key]), str(total_reads_to_CDS), str(length_of_CDS)]
            
            for c in sample_read_codons[site][key].keys():
                codon_c_in_gene = codon_usage[key][c]
                site_with_codon_c = sample_read_codons[site][key][c]
                
                try:
                    normalized_codon_freq = (site_with_codon_c)/((total_reads_to_CDS/length_of_CDS) * codon_c_in_gene)
                    outlist.append(str(normalized_codon_freq))
                except ZeroDivisionError as e:
                    if (site_with_codon_c == total_reads_to_CDS) and (site_with_codon_c == 0):
                        outlist.append(str(0))
                    elif (site_with_codon_c == 0) and (codon_c_in_gene == 0):
                        #outlist.append("K")
                        outlist.append(str(0))
                    elif (site_with_codon_c > 0) and (codon_c_in_gene == 0):
                        #outlist.append("M")
                        outlist.append(str(0))
                    elif (site_with_codon_c > 0) and (total_reads_to_CDS == 0):
                        #outlist.append("B")
                        outlist.append(str(0))
                    
                    else:
                        outlist.append("X")
                        
            
            outfh.write(",".join(outlist) + "\n")
            
            
        outfh2 = open("./output_global_codon_patterns/%s/%s_%s_codon_usage.csv" % (dataset, base_name, site), "w")
        headerlist = ["Gene", "Gene ID", "Length of CDS"]
        for codon in CodonsDict.keys():
            headerlist.append("%s" % codon)
        outfh2.write(",".join(headerlist) + "\n")
        
        for key in sorted(codon_usage.keys()):
            if key[0:5] == "ENSTR" or key[0:4] == "ENSG":
                continue
            outlist = [str(key), str(convert[key]), str(protein_lengths[key])]
            
            for c in CodonsDict.keys():
                outlist.append(str(codon_usage[key][c]))
                
            outfh2.write(",".join(outlist) + "\n")
        
        
    print "Finished with %s at " % base_name, time.ctime()
                
    return averagelines

def calculate_global_codon_patterns_nondivided(dataset, base_name, sample_rpkms, codon_usage, protein_lengths, convert, averagelines):
    sites = ["A", "P", "E", "F"]
    sample_read_codons = {}
    sample_read_totals = {}  # total reads per gene    
    for site in sites:
        sample_read_codons[site] = {}
        
    for site in sites:
        for key in sample_rpkms.keys():
            sample_read_codons[site][key] = CodonsDict.copy()
            sample_read_totals[key] = 0
        
    read_totals = {}
    tot_infh = open("./Annotated_totals_per_read/%s/%s_read_counts.tsv" % (dataset, base_name), "r")
    tot_infh.next()  # skip header
    for line in tot_infh:
        s = line.rstrip().split("\t")
        read_totals[s[0]] = float(s[1])
    tot_infh.close()  
        
        
    bamfile = pysam.AlignmentFile("./Annotated_size_filtered_reads/%s/%s_annotated_sized.bam" % (dataset, base_name), "rb")
    for read in bamfile.fetch(until_eof = True):
        if read.reference_name[0:5] == "ENSTR" or read.reference_name[0:4] == "ENSG":
            continue
        read_gn = read.get_tag("gn")
        read_ac = read.get_tag("ac")  # A-site codon
        read_pc = read.get_tag("pc")  # P-site codon
        read_ec = read.get_tag("ec")  # E-site codon
        try:
            read_fc = read.get_tag("fc")  # "far"-site codon
        except:
            continue
        read_al = read.get_tag("al")  # al for "alignment type", CDS or UTR
        
        if read_al == "CDS":
            sample_read_codons["A"][read_gn][read_ac] += float(1)/read_totals[read.query_name]
            sample_read_codons["P"][read_gn][read_pc] += float(1)/read_totals[read.query_name]
            sample_read_codons["E"][read_gn][read_ec] += float(1)/read_totals[read.query_name]
            try:
                sample_read_codons["F"][read_gn][read_fc] += float(1)/read_totals[read.query_name]
            except:
                continue
            sample_read_totals[read_gn] += float(1)/read_totals[read.query_name]
            
    codon_values = {}
    for site in sites:
        codon_values[site] = {}
    for site in sites:
        for key in CodonsDict.keys():
            codon_values[site][key] = []
            
    for site in sites:
        total_reads_in_high_rpkm_cds = 0
        total_length_high_rpkm_cds = 0
        for key in sample_rpkms.keys():
            if sample_rpkms[key] < 10:
                continue
                
            total_reads_to_CDS = sample_read_totals[key]
            length_of_CDS = protein_lengths[key]    
            total_reads_in_high_rpkm_cds += sample_read_totals[key]
            total_length_high_rpkm_cds += protein_lengths[key]
            
            for c in sample_read_codons[site][key].keys():
                codon_c_in_gene = codon_usage[key][c]
                site_with_codon_c = sample_read_codons[site][key][c]
                
                try:

                    normalized_codon_freq = (site_with_codon_c)/((total_reads_to_CDS/length_of_CDS))
                    codon_values[site][c].append(normalized_codon_freq)
                except ZeroDivisionError as e:
                    if (site_with_codon_c == total_reads_to_CDS) and (site_with_codon_c == 0):
                        codon_values[site][c].append(0)
                    elif (site_with_codon_c == 0) and (codon_c_in_gene == 0):
                        #outlist.append("K")
                        codon_values[site][c].append(0)
                    elif (site_with_codon_c > 0) and (codon_c_in_gene == 0):
                        #outlist.append("M")
                        codon_values[site][c].append(0)
                    elif (site_with_codon_c > 0) and (total_reads_to_CDS == 0):
                        #outlist.append("B")
                        codon_values[site][c].append(0)
                    
                    else:
                        codon_values[site][c].append("X")
        
        
    for site in sites:    
        outfh = open("./output_global_codon_patterns/%s/%s_%s_site_codon_patterns.csv" % (dataset, base_name, site), "w")
        headerlist = ["Gene", "Gene ID", "RPKM", "Total Reads in CDS", "Length of CDS"]
        for codon in CodonsDict.keys():
            headerlist.append("%s" % codon.replace("T","U"))
        outfh.write(",".join(headerlist) + "\n")
        
        averageline = ["All", "All", ">=10", str(total_reads_in_high_rpkm_cds), str(total_length_high_rpkm_cds)]
        for codon in CodonsDict.keys():
            averageline.append(str(sum(codon_values[site][codon])/len(codon_values[site][codon])))
        outfh.write(",".join(averageline) + "\n")
        averagelines[site][base_name] = averageline
        
        for key in sample_rpkms.keys():
            total_reads_to_CDS = sample_read_totals[key]
            length_of_CDS = protein_lengths[key]
            
            outlist = [str(key), str(convert[key]), str(sample_rpkms[key]), str(total_reads_to_CDS), str(length_of_CDS)]
            
            for c in sample_read_codons[site][key].keys():
                codon_c_in_gene = codon_usage[key][c]
                site_with_codon_c = sample_read_codons[site][key][c]
                
                try:
                    #normalized_codon_freq = (site_with_codon_c)/((total_reads_to_CDS/length_of_CDS) * codon_c_in_gene)
                    normalized_codon_freq = (site_with_codon_c)/((total_reads_to_CDS/length_of_CDS))
                    outlist.append(str(normalized_codon_freq))
                except ZeroDivisionError as e:
                    if (site_with_codon_c == total_reads_to_CDS) and (site_with_codon_c == 0):
                        outlist.append(str(0))
                    elif (site_with_codon_c == 0) and (codon_c_in_gene == 0):
                        #outlist.append("K")
                        outlist.append(str(0))
                    elif (site_with_codon_c > 0) and (codon_c_in_gene == 0):
                        #outlist.append("M")
                        outlist.append(str(0))
                    elif (site_with_codon_c > 0) and (total_reads_to_CDS == 0):
                        #outlist.append("B")
                        outlist.append(str(0))
                    
                    else:
                        outlist.append("X")
                        
            
            outfh.write(",".join(outlist) + "\n")
            
            
        outfh2 = open("./output_global_codon_patterns/%s/%s_%s_codon_usage.csv" % (dataset, base_name, site), "w")
        headerlist = ["Gene", "Gene ID", "Length of CDS"]
        for codon in CodonsDict.keys():
            headerlist.append("%s" % codon)
        outfh2.write(",".join(headerlist) + "\n")
        
        for key in sorted(codon_usage.keys()):
            if key[0:5] == "ENSTR" or key[0:4] == "ENSG":
                continue
            outlist = [str(key), str(convert[key]), str(protein_lengths[key])]
            
            for c in CodonsDict.keys():
                outlist.append(str(codon_usage[key][c]))
                
            outfh2.write(",".join(outlist) + "\n")
        
        
    print "Finished with %s at " % base_name, time.ctime()
                
                

    return averagelines
    
    
def write_averages_file(dataset, base_names, averagelines):
    sites = ["A", "P", "E", "F"]

    
    for site in sites:
        outfh_ribo = open("./output_global_codon_patterns/%s/%s_%s_site_ribo_average_codon_patterns_RPKM10.csv" % (dataset, dataset, site), "w")
        outfh_total = open("./output_global_codon_patterns/%s/%s_%s_site_total_average_codon_patterns_RPKM10.csv" % (dataset, dataset, site), "w")
        
        headerlist = ["Dataset", "RPKM Cutoff", "Total Reads included", "Length of CDSs included"]
        for codon in CodonsDict.keys():
            headerlist.append("%s" % codon.replace("T","U"))
                
        outfh_ribo.write(",".join(headerlist) + "\n")
        outfh_total.write(",".join(headerlist) + "\n")
        for base_name in base_names:
            if base_name[1] == "R":
                #print repr(averagelines[site][base_name])
                outfh_ribo.write(base_name + "," + ",".join(averagelines[site][base_name][2:]) + "\n")
                
            if base_name[1] == "T":
                #print repr(averagelines[site][base_name])
                outfh_total.write(base_name + "," + ",".join(averagelines[site][base_name][2:]) + "\n")
                
                
        outfh_ribo.close()
        outfh_total.close()
        
        
def write_averages_file_nondivided(dataset, base_names, averagelines):
    sites = ["A", "P", "E", "F"]

    
    for site in sites:
        outfh_ribo = open("./output_global_codon_patterns_nondivided/%s/%s_%s_site_ribo_average_codon_patterns_nondivided_RPKM10.csv" % (dataset, dataset, site), "w")
        outfh_total = open("./output_global_codon_patterns_nondivided/%s/%s_%s_site_total_average_codon_patterns_nondivided_RPKM10.csv" % (dataset, dataset, site), "w")
        
        headerlist = ["Dataset", "RPKM Cutoff", "Total Reads included", "Length of CDSs included"]
        for codon in CodonsDict.keys():
            headerlist.append("%s" % codon.replace("T","U"))
                
        outfh_ribo.write(",".join(headerlist) + "\n")
        outfh_total.write(",".join(headerlist) + "\n")
        for base_name in base_names:
            if base_name[1] == "R":
                outfh_ribo.write(base_name + "," + ",".join(averagelines[site][base_name][2:]) + "\n")
                
            if base_name[1] == "T":
                outfh_total.write(base_name + "," + ",".join(averagelines[site][base_name][2:]) + "\n")
                
                
        outfh_ribo.close()
        outfh_total.close()    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Filter the annotated reads to the appropriate size for the dataset.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the Annotated_reads folder (S1, S2, etc)")
    parser.add_argument("protein_filepath", help = "The path to the protein file (e.g. BowtieIndex/S12_all_protein.fa) (needed because it compares CDS translations to make sure it gets the frame right")
    parser.add_argument("cds_filepath", help = "The path to the CDS file (e.g. BowtieIndex/S12_all_CDS.fa)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    protein_sequences, protein_lengths = read_proteins(args.protein_filepath)
    codon_usage = calculate_codon_usage(args.cds_filepath, protein_sequences)
    coordinates_dict = get_start_end_coords(args.dataset)
    convert = ensembl_ID_converter(args.dataset)
    sites = ["A", "P", "E", "F"]
        
    data_files = sorted([x for x in os.listdir("./Annotated_size_filtered_reads/%s/" % args.dataset) if x[-4:] == ".bam"])
    base_names = []
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[0:-2])
        base_names.append(base_name)
        
    print base_names
        
    averagelines = {}
    for site in sites:
        averagelines[site] = {}
        
    averagelines_nondivided = {}
    for site in sites:
        averagelines_nondivided[site] = {}        
        
    for f in data_files:
        base_name = "_".join(f.split(".")[0].split("_")[0:-2])
        sample_rpkms = read_rpkm_values(args.dataset, base_name, 0)  # can use a cutoff other than 0 if desired
        averagelines = calculate_global_codon_patterns(args.dataset, base_name, sample_rpkms, codon_usage, protein_lengths, convert, averagelines)
        averagelines_nondivided = calculate_global_codon_patterns_nondivided(args.dataset, base_name, sample_rpkms, codon_usage, protein_lengths, convert, averagelines_nondivided)
        
    write_averages_file(args.dataset, base_names, averagelines)
    write_averages_file_nondivided(args.dataset, base_names, averagelines_nondivided)
    

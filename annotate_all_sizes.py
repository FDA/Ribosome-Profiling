#!/usr/bin/env python  

import sys
import os
import math
import subprocess
import argparse
import time
from collections import OrderedDict
import pysam


def create_output_folders(dataset):
    if not os.path.exists("./Annotated_reads"):
        os.mkdir("./Annotated_reads")
    if not os.path.exists("./Annotated_reads/%s" % args.dataset):
        os.mkdir("./Annotated_reads/%s" % args.dataset)
        
class DataSetOffsetInfo(object):
    def __init__(self, offset):
        self.offset = offset
        
    def add_offset(self, key, value):
        self.offset[key] = value
        
def define_dataset_offset_info():
    global offsetinfo
    offsetinfo = {}
    offsetinfo["S12"] = DataSetOffsetInfo({20: 15, 21: 15, 22: 15, 23: 15, 24: 15, 25: 15, 26: 15, 27: 15, 28: 15, 29: 15, 30: 15, 31: 15, 32: 15})
    for i in range(33, 126):
        offsetinfo["S12"].add_offset(i, 15)
        
class Transcript(object):
    def __init__(self, ident, start, end):
        self.transcript_id = ident
        self.start = start
        self.end = end
        
def transcript_start_end(dataset):
    transcripts = {}
    infh = open("%s_transcript_CDS_start_end_locations.txt" % dataset, "r")

    # skip header line
    next(infh)

    for line in infh:
        s = line.rstrip().split()
        tid, start, end = s[0], int(s[1])-1, int(s[2])-1
        transcripts[tid] = Transcript(tid, start, end)
        
    return transcripts
        
        
def annotate_reads(dataset, base_name, transcripts):
    
     
    bamfile = pysam.AlignmentFile("./Hisat_output/%s/%s/%s_sorted.bam" % (dataset, base_name, base_name), "rb")

    outsam = pysam.AlignmentFile("./Annotated_reads/%s/annotated_%s.bam" % (dataset, base_name), mode = "wb", template = bamfile)
     
    for read in bamfile.fetch():
        gene_name = read.reference_name
        
        if (read.query_length not in offsetinfo[dataset].offset.keys()):
            continue
        
        gene_length = transcripts[gene_name].end - transcripts[gene_name].start +1
        offset_value = offsetinfo[dataset].offset[read.query_length]
        
        if ((read.reference_start+offset_value) >= transcripts[gene_name].start):
            if (read.reference_start+offset_value+2) >= transcripts[gene_name].end:
                # if the A-site is past the CDS
                align_type = "UTR"
                align_type_more = "3UTR"
            else:
                # if the A-site is in the CDS
                align_type = "CDS" 
                align_type_more = "CDS"
                
        else:
            # if the A-site is before the CDS
            align_type = "UTR"
            align_type_more = "5UTR"

        if align_type == "CDS":
            if (read.reference_start) >= transcripts[gene_name].start:
                frame = (read.reference_start - transcripts[gene_name].start) % 3
                pairs = read.get_aligned_pairs(matches_only = True, with_seq = True)
                if  len(pairs) < 20:
                    continue
                if frame == 0:                    
                    # A-site codon
                    codon = (pairs[offset_value][2] + pairs[offset_value+1][2] + pairs[offset_value+2][2]).upper()

                    # P-site codon
                    pcodon = (pairs[offset_value-3][2] + pairs[offset_value-2][2] + pairs[offset_value-1][2]).upper()

                    # E-site codon
                    ecodon = (pairs[offset_value-6][2] + pairs[offset_value-5][2] + pairs[offset_value-4][2]).upper()
                    try: 
                        # "Far" codon, on opposite side of A from P
                        fcodon = (pairs[offset_value+3][2] + pairs[offset_value+4][2] + pairs[offset_value+5][2]).upper()
                    except: 
                        pass
                    position = pairs[offset_value][1]
                    pposition = pairs[offset_value-3][1]
                    eposition = pairs[offset_value-6][1]
                    try: 
                        fposition = pairs[offset_value+3][1]
                    except: 
                        pass

                elif frame == 1:
                    codon = (pairs[offset_value-1][2] + pairs[offset_value][2] + pairs[offset_value+1][2]).upper()
                    pcodon = (pairs[offset_value-4][2] + pairs[offset_value-3][2] + pairs[offset_value-2][2]).upper()
                    ecodon = (pairs[offset_value-7][2] + pairs[offset_value-6][2] + pairs[offset_value-5][2]).upper()
                    try: 
                        fcodon = (pairs[offset_value+2][2] + pairs[offset_value+3][2] + pairs[offset_value+4][2]).upper()
                    except: 
                        pass
                    position = pairs[offset_value-1][1]
                    pposition = pairs[offset_value-4][1]
                    eposition = pairs[offset_value-7][1]
                    try: 
                        fposition = pairs[offset_value+2][1]
                    except: 
                        pass
                elif frame == 2:
                    codon = (pairs[offset_value+1][2] + pairs[offset_value+2][2] + pairs[offset_value+3][2]).upper()
                    pcodon = (pairs[offset_value-2][2] + pairs[offset_value-1][2] + pairs[offset_value][2]).upper()
                    ecodon = (pairs[offset_value-5][2] + pairs[offset_value-4][2] + pairs[offset_value-3][2]).upper()
                    try: 
                        fcodon = (pairs[offset_value+4][2] + pairs[offset_value+5][2] + pairs[offset_value+6][2]).upper()
                    except: 
                        pass
                    position = pairs[offset_value+1][1]
                    pposition = pairs[offset_value-2][1]
                    eposition = pairs[offset_value-5][1]
                    try: 
                        fposition = pairs[offset_value+4][1]
                    except: 
                        pass
        
            
                local_position = position - transcripts[gene_name].start
                plocal_position = pposition - transcripts[gene_name].start
                elocal_position = eposition - transcripts[gene_name].start
                try: 
                    flocal_position = fposition - transcripts[gene_name].start
                except: 
                    pass
                AA_position = int(math.floor(float(local_position)/3))
                
            elif (read.reference_start) < transcripts[gene_name].start:
                offset_from_start_codon = transcripts[gene_name].start - read.reference_start
                frame_wrong_direction = offset_from_start_codon%3
                if frame_wrong_direction != 0:
                    frame = 3 - frame_wrong_direction
                else:
                    frame = 0

                pairs = read.get_aligned_pairs(matches_only = True, with_seq = True)
                if  len(pairs) < 20:
                    continue
                if frame == 0:
                    # A-site codon
                    codon = (pairs[offset_value][2] + pairs[offset_value+1][2] + pairs[offset_value+2][2]).upper()

                    # P-site codon
                    pcodon = (pairs[offset_value-3][2] + pairs[offset_value-2][2] + pairs[offset_value-1][2]).upper()

                    # E-site codon
                    ecodon = (pairs[offset_value-6][2] + pairs[offset_value-5][2] + pairs[offset_value-4][2]).upper()
                    try: 
                        # "Far" codon, on opposite side of A from P
                        fcodon = (pairs[offset_value+3][2] + pairs[offset_value+4][2] + pairs[offset_value+5][2]).upper()
                    except: 
                        pass
                    position = pairs[offset_value][1]
                    pposition = pairs[offset_value-3][1]
                    eposition = pairs[offset_value-6][1]
                    try: 
                        fposition = pairs[offset_value+3][1]
                    except: 
                        pass
                elif frame == 1:
                    codon = (pairs[offset_value-1][2] + pairs[offset_value][2] + pairs[offset_value+1][2]).upper()
                    pcodon = (pairs[offset_value-4][2] + pairs[offset_value-3][2] + pairs[offset_value-2][2]).upper()
                    ecodon = (pairs[offset_value-7][2] + pairs[offset_value-6][2] + pairs[offset_value-5][2]).upper()
                    try: 
                        fcodon = (pairs[offset_value+2][2] + pairs[offset_value+3][2] + pairs[offset_value+4][2]).upper()
                    except: 
                        pass
                    position = pairs[offset_value-1][1]
                    pposition = pairs[offset_value-4][1]
                    eposition = pairs[offset_value-7][1]
                    try: 
                        fposition = pairs[offset_value+2][1]
                    except: 
                        pass
                elif frame == 2:
                    codon = (pairs[offset_value+1][2] + pairs[offset_value+2][2] + pairs[offset_value+3][2]).upper()
                    pcodon = (pairs[offset_value-2][2] + pairs[offset_value-1][2] + pairs[offset_value][2]).upper()
                    ecodon = (pairs[offset_value-5][2] + pairs[offset_value-4][2] + pairs[offset_value-3][2]).upper()
                    try: 
                        fcodon = (pairs[offset_value+4][2] + pairs[offset_value+5][2] + pairs[offset_value+6][2]).upper()
                    except: 
                        pass
                    position = pairs[offset_value+1][1]
                    pposition = pairs[offset_value-2][1]
                    eposition = pairs[offset_value-5][1]
                    try: 
                        fposition = pairs[offset_value+4][1]
                    except: 
                        pass

                local_position = position - transcripts[gene_name].start
                plocal_position = pposition - transcripts[gene_name].start
                elocal_position = eposition - transcripts[gene_name].start
                try: 
                    flocal_position = fposition - transcripts[gene_name].start
                except: 
                    pass
                AA_position = int(math.floor(float(local_position)/3))
            
            
        else:
            codon = "NA"
            pcodon = "NA"
            ecodon = "NA"
            try: 
                fcodon = "NA"
            except: 
                pass
            local_position = -1
            plocal_position = -1
            elocal_position = -1
            try: 
                flocal_position = -1
            except: 
                pass
            AA_position = -1
            frame = -1
            
        read.set_tag("gn", gene_name, "Z")
        read.set_tag("gl", gene_length, "i")
        read.set_tag("ac", codon, "Z")
        read.set_tag("pc", pcodon, "Z")
        read.set_tag("ec", ecodon, "Z")
        try: 
            read.set_tag("fc", fcodon, "Z")
        except: 
            pass
        read.set_tag("ov", offset_value, "i")
        read.set_tag("lp", local_position, "i")
        read.set_tag("as", local_position, "i")
        read.set_tag("ps", plocal_position, "i")
        read.set_tag("es", elocal_position, "i")
        try: 
            read.set_tag("fs", flocal_position, "i")
        except: 
            pass
        read.set_tag("al", align_type, "Z")
        read.set_tag("ut", align_type_more, "Z")
        read.set_tag("pp", AA_position, "i")
        read.set_tag("rf", frame, "i")
        
        outsam.write(read)
            
    print("Finished with %s at" % base_name, time.ctime())
        
        
def index_new_bamfile(dataset, base_name):   
    pysam.index("./Annotated_reads/%s/annotated_%s.bam" % (dataset, base_name))
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Annotate features of reads in a dataset after adapter trimming, rRNA+tRNA decontamination, and alignment.")
    parser.add_argument("dataset", help = "The folder you want to look in for the data, inside the HisatOutput folder (S1, S2, etc)")
    args = parser.parse_args()
    
    create_output_folders(args.dataset)
    data_folders = sorted(os.listdir("./Hisat_output/%s/" % args.dataset))
    define_dataset_offset_info()
    transcripts = transcript_start_end(args.dataset)
    
    for f in data_folders:
        base_name = f
        annotate_reads(args.dataset, base_name, transcripts)
        index_new_bamfile(args.dataset, base_name)

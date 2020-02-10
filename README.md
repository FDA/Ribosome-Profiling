The following are instructions for setting up and running the ribosome profiling
analysis pipeline as described in: 'Effects of codon optimization on
coagulation factor IX translation and structure: Implications for protein and
gene therapies' (Alexaki et al., 2019)

Raw data for the above paper are accessible from NCBI SRA Database
https://www.ncbi.nlm.nih.gov/sra/
Project accession: PRJNA591214

They may be downloaded using the NCBI SRA toolkit
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

The following are the specific run accession numbers and original filenames for
the datafiles used in the above paper:

SRR10513636 - 1R_CAGATC_L002_R1_001.fastq.gz
SRR10513635 - 1T_ATCACG_L002_R1_001.fastq.gz
SRR10513632 - 2R_ACTTGA_L002_R1_001.fastq.gz
SRR10513631 - 2T_CGATGT_L002_R1_001.fastq.gz
SRR10513630 - 3R_GATCAG_L002_R1_001.fastq.gz
SRR10513629 - 3T_TTAGGC_L002_R1_001.fastq.gz
SRR10513628 - 4R_TAGCTT_L002_R1_001.fastq.gz
SRR10513627 - 4T_TGACCA_L002_R1_001.fastq.gz
SRR10513626 - 5R_GGCTAC_L002_R1_001.fastq.gz
SRR10513625	- 5T_ACAGTG_L002_R1_001.fastq.gz
SRR10513634 - 6R_CTTGTA_L002_R1_001.fastq.gz
SRR10513633 - 6T_GCCAAT_L002_R1_001.fastq.gz

********************************************************************************

If running this pipeline on other datasets, it is advisable that the user check
the data for contamination with common bacterial or viral sequences, depending
on how the data were gathered. Additional decontamination steps may need to be
taken at the user's discretion.

********************************************************************************

All scripts were written by John Athey while working in Dr. Chava
Kimchi-Sarfaty's laboratory at the FDA, White Oak, MD, USA.

This pipeline represents Version 1.2

Prerequisites:

Python 2.7 (https://www.python.org/)
GFF Utilities (http://ccb.jhu.edu/software/stringtie/gff.shtml)
Bowtie (http://bowtie-bio.sourceforge.net/index.shtml)
TopHat (https://ccb.jhu.edu/software/tophat/index.shtml)
FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
Samtools (http://www.htslib.org/)

Folder structure should be setup as follows:

./Ribosome_profiling/
    /Raw_data/
    /BowtieIndex/

The following files should be placed in ./Ribosome_profiling/BowtieIndex/

    annotation_filter.py
    contaminant_sequences.fa
    ensembl_CDS_lengths.py
    ensembl_to_gene_id.py
    ensembl_transcript_lengths.py
    F9_opt1_construct_100bpUTRs.fasta
    F9_WT_construct_100bpUTRs.fasta
    opt1F9c_annotation.gff3
    WTF9c_annotation.gff3
    
All other Python scripts should be placed in ./Ribosome_profiling/

##### Instructions for building custom Bowtie index: #####
##### Must be completed before running main pipeline #####

Download GRCh37.p13.genome.fa and gencode.v19.annotation.gff3 to 
'./Ribosome_profiling/BowtieIndex/' GRCh37.p13.genome.fa and 
gencode.v19.annotation.gff3 were used in: 'Effects of codon optimization on
coagulation factor IX translation and structure: Implications for protein and
gene therapies' (Alexaki et al., 2019). Alternatively, these steps could also be
followed using the newest (as of November 20, 2019) human assembly and
annotations, ie, GRCh38.p13.genome.fa and gencode.v32.annotation.gff3 
(see https://www.gencodegenes.org/human/). If using a different assembly and 
annotations, adjust accordingly in the script 'annotation_filter.py'. The 
current version of the script has 'gencode.v19.annotation.gff3' hardcoded as the
input and 'filtered_gencode.v19.annotation.gff3' hardcoded as the output. Also, 
adjust any of the commands for setting up bowtie index that make reference to 
the assembly or annotation version. 'contaminant_sequences.fa' will also need to 
be replaced with a version that corresponds to the assembly and annotations 
being used.

NOTE: THIS IS NOT A COMMAND. S12 is the dataset name we have used, but it can be 
changed. If changed, make sure to adjust the argument given to subsequent 
scripts in the pipeline. The dataset argument will tell the script where to look
for the input files. Pipeline begins with the raw data files (*.gz) in 
./Ribosome_profiling/Raw_data/S12/

# denotes a comment, not to be run as a command #

The following commands can be run in a unix/linux terminal. Portions in '()' are
not part of command.

If the following commands are run verbatim, they should be executed from 
'./Ribosome_profiling/BowtieIndex/'
To run without any alterations, place raw data files in 
'./Ribosome_profiling/Raw_data/S12/'

##### BOWTIE INDEX SETUP COMMANDS #####

python annotation_filter.py (-> filtered_gencode.v19.annotation.gff3)
grep -v "gene_name=F9" filtered_gencode.v19.annotation.gff3 > annotation_without_WT_F9.gff3     #Remove WT F9 from annotations
cp GRCh37.p13.genome.fa genome_plus_F9_constructs.fa #Creates new genome that constructs will be added to
cat F9_WT_construct_100bpUTRs.fasta >> genome_plus_F9_constructs.fa #Add constructs to genome
cat F9_opt1_construct_100bpUTRs.fasta >> genome_plus_F9_constructs.fa
grep -v "gene_name=F9" filtered_gencode.v19.annotation.gff3 > annotation_mWTF9_pWTF9c.gff3 #All annotations minus WTF9
cat WTF9c_annotation.gff3 >> annotation_mWTF9_pWTF9c.gff3 #Add WTF9 construct annotation
grep -v "gene_name=F9" filtered_gencode.v19.annotation.gff3 > annotation_mWTF9_popt1F9c.gff3 #All annotations minus WTF9
cat opt1F9c_annotation.gff3 >> annotation_mWTF9_popt1F9c.gff3 #Add optimized F9 construct annotation
gffread -W -w transcripts_mWTF9_pWTF9c.fa -x CDS_mWTF9_pWTF9c.fa -y proteins_mWTF9_pWTF9c.fa -C -g genome_plus_F9_constructs.fa annotation_mWTF9_pWTF9c.gff3
gffread -W -w transcripts_mWTF9_popt1F9c.fa -x CDS_mWTF9_popt1F9c.fa -y proteins_mWTF9_popt1F9c.fa -C -g genome_plus_F9_constructs.fa annotation_mWTF9_popt1F9c.gff3
bowtie-build transcripts_mWTF9_pWTF9c.fa transcripts_mWTF9_pWTF9c  #Build Reference
bowtie-build transcripts_mWTF9_popt1F9c.fa transcripts_mWTF9_popt1F9c
bowtie-build contaminant_sequences.fa contaminant_sequences
cat annotation_without_WT_F9.gff3 WTF9c_annotation.gff3 opt1F9c_annotation.gff3 > S12_all_annotations.gff3
gffread -W -w S12_all_transcripts.fa -x S12_all_CDS.fa -y S12_all_protein.fa -C -g genome_plus_F9_constructs.fa S12_all_annotations.gff3
python ensembl_to_gene_id.py S12
python ensembl_transcript_lengths.py S12
python ensembl_CDS_lengths.py S12
cp S12* ..  (to ~/Ribosome_profiling)

NOTE: THIS IS NOT A COMMAND. S12 is the dataset name we have used, but it can be
changed. If changed, make sure to adjust the argument given to subsequent
scripts in the pipeline. Pipeline begins with the raw data files (*.gz) in
./Ribosome_profiling/Raw_data/S12/

******************************************************
*                                                    *
*    MAIN PIPELINE RUN FROM ./Ribosome_profiling/    *
*                                                    *
******************************************************

    #Lines that start with '#' are instructions to the user

python trim_adapters.py S12
python decontaminate_reads.py S12
python calculate_length_distribution_pre_alignment.py S12
    #If running for S12, no further action is required. Otherwise, open 
    #'initial_read_alignment.py' and adjust 'define_dataset_info()' accordingly.
python initial_read_alignment.py S12
python samtools_sort_index.py S12
python util_parse_CDS_start_end.py S12
python calculate_length_distribution_post_alignment.py S12
python calculate_p_site_offset.py S12
python calcualte_stop_codon_p_site_offet.py S12
    #If running for S12, no further action is required. Otherwise, open 
    #'annotate_all_sizes.py' and adjust offset info for your dataset based on 
    #output from previous two scripts.
python annotate_all_sizes.py S12
python util_calculate_annotated_totatls_per_read.py S12
python calculate_frame_distribution.py S12
python calculate_UTR_CDS_distribution.py S12
    #Open 'annotated_size_filter.py' and decided lengths that 
    #represent true RPFs based on validation output. 
    #Adjust 'define_dataset_sizes()' accordingly.
python annotated_size_filter.py S12
python calculate_rpkm.py S12
python calculate_global_codon_patterns.py S12 S12_all_protein.fa S12_all_CDS.fa
    #Open 'calculate_gene_of_interest.py' and ensure sizes are defined 
    #correctly. The following script is run with gene of interest as an argument. 
python calculate_gene_of_interest.py S12 F9WT       
    #For other genes, replace F9WT. 
    #(F9opt1, ACTB, GAPDH, or any other gene of interest)
    #Open 'calculate_mRNA_coverage.py' and ensure 'special_cases' are correct. 
    #For S12 they are 'F9WT' and 'F9opt1'.
python calculate_mRNA_coverage.py S12 F9WT #Same comment as above
python util_gene_of_interest_formatter.py S12 F9WT #Same comment as above

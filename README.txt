The following are instruction for setting up and running the ribosome profiling
analysis pipeline as described in: 'Effects of codon optimization on 
coagulation factor IX translation and structure: Implications for protein and
gene therapies' Alexaki et al. 2019, with modifications. In the aforementioned
manuscript, TopHat (version 2.0.9) was used in the pipeline for the alignment
step. In the present pipeline, HISAT2 (version 2.1.0) is used for the alignment
step.

All scripts were written by John Athey while working in Dr. Chava
Kimchi-Sarfaty's laboratory at the FDA, White Oak, MD, USA.

This pipeline represents Version 2.2

The following prerequisites (version tested) must be met by the user before
executing the pipeline:

Python 3.7 (3.7.6) (https://www.python.org/)
    libraries:
    pysam (0.15.3) (https://github.com/pysam-developers/pysam)
    biopython (1.77) (https://biopython.org/)
GFF Utilities (gffread v0.12.1) (http://ccb.jhu.edu/software/stringtie/gff.shtml)
Bowtie (1.0.0) (http://bowtie-bio.sourceforge.net/index.shtml)
HISAT2 (2.1.0) (https://ccb.jhu.edu/software/hisat2/manual.shtml)
FASTX-Toolkit (0.0.14) (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
Samtools (1.7 using htslib 1.7) (http://www.htslib.org/)


Instructions for building custom HISAT index:

From https://www.gencodegenes.org/human/release_19.html download and unzip the
human genome assembly on all regions (GRCh37.p13.genome.fa.gz) and
comprehensive gene annotations (gencode.v19.annotation.gff3.gz) to 
'./Ribosome_profiling/HisatIndex/'. Please note that GRCh37.p13.genome.fa and
gencode.v19.annotation.gff3 were used in: 'Effects of codon optimization on
coagulation factor IX translation and structure: Implications for protein and
gene therapies' Alexaki et al. 2019. Alternatively, these steps could also be
followed using the newest (as of August 3, 2020) human genome assembly on all
regions and comprehensive gene annotations on the reference chromosomes only,
ie, GRCh38.p13.genome.fa and gencode.v34.annotation.gff3
(see https://www.gencodegenes.org/human/). If using a different assembly and
annotations, adjust accordingly in the script 'annotation_filter.py'
(see comments within script). The current version of the script has
'gencode.v19.annotation.gff3' hardcoded as the input and
'filtered_gencode.v19.annotation.gff3' hardcoded as the output. Also, adjust
any of the commands for setting up Hisat index that make reference to the
assembly or annotation version (see comments in 'build_hisat_index.sh'.
'contaminant_sequences.fa' will also need to be replaced with a version that
corresponds to the assembly and annotations being used. In the current version,
the dataset is defined in 'build_hisat_index.sh' as 'S12'. If it is changed,
the './Ribosome_profiling/Raw_data/S12/' folder mentioned below must be changed
to match the dataset name as well.

To build the HISAT index, run the bash script 'build_hisat_index.sh' from
'./Ribosome_profiling/':

bash build_hisat_index.sh



Download .fastq.gz files from bioproject PRJNA591214 (12 in total)
(https://www.ncbi.nlm.nih.gov/bioproject/591214) to
'./Ribosome_profiling/Raw_data/S12/'. As mentioned above, if the user chooses
a different dataset name, please replace 'S12' with the dataset name, and
adjust the dataset name defined in 'RP_analysis_pipeline.sh'.

Run the bash script 'RP_analysis_pipeline.sh' from './Ribosome_profiling/':

bash RP_analysis_pipeline.sh

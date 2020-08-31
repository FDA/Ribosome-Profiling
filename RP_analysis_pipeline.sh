#directory paths
cwd=$(pwd)
dataset="S12"

python trim_adapters.py $dataset
python decontaminate_reads.py $dataset
python calculate_length_distribution_pre_alignment.py $dataset

###    If running for S12, no further action is required. Otherwise, open 'initial_read_alignment.py' and adjust 'define_dataset_info()' accordingly.

python initial_read_alignment.py $dataset
python samtools_sort_index.py $dataset
python util_parse_CDS_start_end.py $dataset
python calculate_length_distribution_post_alignment.py $dataset
python calculate_p_site_offset.py $dataset
python calculate_stop_codon_p_site_offset.py $dataset

###    If running for S12, no further action is required. Otherwise, open 'annotate_all_sizes.py' and adjust offset info for your dataset based on output from previous two scripts.

python annotate_all_sizes.py $dataset
python util_calculate_annotated_totals_per_read.py $dataset
python calculate_frame_distribution.py $dataset
python calculate_UTR_CDS_distribution.py $dataset

###    Open 'annotated_size_filter.py' and decided lengths that represent true RPFs based on validation output. Adjust 'define_dataset_sizes()' accordingly.

python annotated_size_filter.py $dataset
python calculate_rpkm.py $dataset

###    Open 'calculate_gene_of_interest.py' and ensure sizes are defined correctly. The following script is run with gene of interest as an argument. 

###    Open 'calculate_mRNA_coverage.py' and ensure 'special_cases' are correct. For S12 they are 'F9WT' and 'F9opt1'.

### For other genes of interest, replace below
for gene in F9WT F9opt1 ACTB GAPDH
do
  python calculate_gene_of_interest.py $dataset $gene
  python calculate_mRNA_coverage.py $dataset $gene
  python util_gene_of_interest_formatter.py $dataset $gene
done



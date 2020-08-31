#directory paths
cwd=$(pwd)
hisatd="${cwd}/HisatIndex"
dataset="S12"


cd ${hisatd}

python annotation_filter.py
#The following command can be adjusted if using a different gencode version
grep -v "gene_name=F9" filtered_gencode.v19.annotation.gff3 > annotation_without_WT_F9.gff3

#The following command can be adjusted if using a different assembly
cp GRCh37.p13.genome.fa genome_plus_F9_constructs.fa

cat F9_WT_construct_100bpUTRs.fasta >> genome_plus_F9_constructs.fa
cat F9_opt1_construct_100bpUTRs.fasta >> genome_plus_F9_constructs.fa

#The following grep commands can be adjusted if using a different gencode verison
grep -v "gene_name=F9" filtered_gencode.v19.annotation.gff3 > annotation_mWTF9_pWTF9c.gff3
cat WTF9c_annotation.gff3 >> annotation_mWTF9_pWTF9c.gff3
grep -v "gene_name=F9" filtered_gencode.v19.annotation.gff3 > annotation_mWTF9_popt1F9c.gff3
cat opt1F9c_annotation.gff3 >> annotation_mWTF9_popt1F9c.gff3

gffread -W -w transcripts_mWTF9_pWTF9c.fa -x CDS_mWTF9_pWTF9c.fa -y proteins_mWTF9_pWTF9c.fa -C -g genome_plus_F9_constructs.fa annotation_mWTF9_pWTF9c.gff3
gffread -W -w transcripts_mWTF9_popt1F9c.fa -x CDS_mWTF9_popt1F9c.fa -y proteins_mWTF9_popt1F9c.fa -C -g genome_plus_F9_constructs.fa annotation_mWTF9_popt1F9c.gff3
hisat2-build transcripts_mWTF9_pWTF9c.fa transcripts_mWTF9_pWTF9c
hisat2-build transcripts_mWTF9_popt1F9c.fa transcripts_mWTF9_popt1F9c
bowtie-build contaminant_sequences.fa contaminant_sequences
cat annotation_without_WT_F9.gff3 WTF9c_annotation.gff3 opt1F9c_annotation.gff3 > ${dataset}_all_annotations.gff3
gffread -W -w ${dataset}_all_transcripts.fa -x ${dataset}_all_CDS.fa -y ${dataset}_all_protein.fa -C -g genome_plus_F9_constructs.fa ${dataset}_all_annotations.gff3
python ensembl_to_gene_id.py ${dataset}
python ensembl_transcript_lengths.py ${dataset}
python ensembl_CDS_lengths.py ${dataset}
cp ${dataset}* ..

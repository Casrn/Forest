#!/bin/bash
#title: ITS processing with Qiime2(DADA2)
#author: Yaping Liu
#The sequence data is a set of Novaseq-sequenced paired-end fastq files. 
#primer: ITS2: ITS3-ITS4 GCATCGATGAAGAACGCAGC-TCCTCCGCTTATTGATATGC

#---
#last updated: 2024-05-09
#---

# Before using Qiime2, rename the sequence file names
## 1-60：a1-a60
#for i in {1..60}
#do 
	#mv "${i}_R1.fastq.gz" "a${i}_R1.fastq.gz"
	#mv "${i}_R2.fastq.gz" "a${i}_R2.fastq.gz"
#done 

## 61-120：b-3-b57(b-3, b-2, b-1, b1,b2,....b57)
#for i in {61..120}
#do 
	#mv "${i}_R1.fastq.gz" "b$((i-63))_R1.fastq.gz"
	#mv "${i}_R2.fastq.gz" "b$((i-63))_R2.fastq.gz"
#done

## 121-180：c1-c60
#for i in {121..180}
#do 
	#mv "${i}_R1.fastq.gz" "c$((i-120))_R1.fastq.gz"
	#mv "${i}_R2.fastq.gz" "c$((i-120))_R2.fastq.gz"
#done

# Activate Qiime2 environment
source /usr/local/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2023.2

# Set path and sequence of primer
SEQ="/home/dell/02Data/Forest/ITS_20240604"
OUT="/home/liuyaping/GW/northforest_rawdata/ITS"
DATABASE="/home/dell/00Localdb"
# Create Manifest file
cd $SEQ
ls *fq.gz | \
    awk -F'_null' -v P="$(pwd)" 'BEGIN { print "sample-id,absolute-filepath,direction" } \
    {if(NR%2==0){print $1","P"/"$0",reverse"}else{print $1","P"/"$0",forward"}}' > $OUT/sample_dataset.csv
cd $OUT

# using fastqc to check the quality of sequences 
#mkdir $OUT/fastqc_out
#fastqc -t 4 /home/dell/02Data/MuUs/MuUs_180sample_202309/ITS/* \
#-o /home/liuyaping/MuUs/data_processing/qiime2/fastqc_out_ITS

# Import the data #3 mins
time qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path $OUT/sample_dataset.csv \
	--input-format PairedEndFastqManifestPhred33  \
	--output-path $OUT/paired-end-demux.qza

# Trim primer and adapter
# 515F-806R
time qiime cutadapt trim-paired \
	--p-cores 30 \
	--p-discard-untrimmed \
	--p-error-rate 0.1 \
	--i-demultiplexed-sequences $OUT/paired-end-demux.qza\
	--p-front-f AACTTTYRRCAAYGGATCWCT \
	--p-front-r AGCCTCCGCTTATTGATATGCTTAART \
	--o-trimmed-sequences $OUT/paired-end-demux-trimmed.qza\
	--verbose

# Visualize quality of reads
# https://view.qiime2.org
## before trimming
time qiime demux summarize \
	--i-data $OUT/paired-end-demux.qza \
	--o-visualization $OUT/paired-end-demux.qzv
## after trimming
time qiime demux summarize \
	--i-data $OUT/paired-end-demux-trimmed.qza \
	--o-visualization $OUT/paired-end-demux-trimmed.qzv

# Denoise with DADA2
mkdir $OUT/dada2
time qiime dada2 denoise-paired \
	--i-demultiplexed-seqs $OUT/paired-end-demux-trimmed.qza \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-trunc-len-f 180 \
	--p-trunc-len-r 200 \
	--p-max-ee-f 2 \
	--p-max-ee-r 2 \
	--p-chimera-method consensus \
	--o-denoising-stats $OUT/dada2/stats-dada2.qza \
	--o-table $OUT/dada2/fungi_count.qza \
	--o-representative-sequences $OUT/dada2/fungi_rep_seqs.qza \
	--verbose

# Visualize denoized results
# https://view.qiime2.org
## table
time qiime feature-table summarize \
	--i-table $OUT/dada2/fungi_count.qza\
	--o-visualization $OUT/dada2/fungi_count.qzv
## rep-seq
time qiime feature-table tabulate-seqs \
	--i-data $OUT/dada2/fungi_rep_seqs.qza\
	--o-visualization $OUT/dada2/fungi_rep_seqs.qzv

# Assign taxonomy
## self-trained Silva Naive Bayes classifier
mkdir $OUT/classified_sequences
time qiime feature-classifier classify-sklearn \
	--i-classifier $DATABASE/unite-ver9-dynamic-classifier-03.16.2023.qza \
	--i-reads $OUT/dada2/fungi_rep_seqs.qza \
	--o-classification $OUT/classified_sequences/classification.qza

# Visualize taxonomic classification
# https://view.qiime2.org
## Visualize taxonomic classification
time qiime metadata tabulate \
	--m-input-file $OUT/classified_sequences/classification.qza \
	--o-visualization $OUT/classified_sequences/classification.qzv
## Visualize feature stats table
time qiime metadata tabulate \
	--m-input-file $OUT/dada2/fungi_count.qza \
	--o-visualization $OUT/dada2/stats-dada2.qzv
	
# Export the data 
## taxonomic table 
qiime tools export \
	--input-path $OUT/classified_sequences/classification.qza \
	--output-path $OUT/final_export
mv $OUT/final_export/taxonomy.tsv $OUT/final_export/fungi_taxonomy.tsv
## feature table 
qiime tools export \
	--input-path $OUT/dada2/fungi_count.qza \
	--output-path $OUT/final_export
biom convert \
	-i $OUT/final_export/feature-table.biom \
    -o $OUT/final_export/fungi_count.txt \
    --to-tsv	
# Build phylogeny
## qiime fragment-insertion tree with a reference database
#qiime fragment-insertion sepp \
	#--i-representative-sequences $OUT/classified_sequences/classification.qza \
	#--i-reference-database $DATABASE/sepp-refs-silva-128.qza \
	#--o-tree $OUT/insertion-tree.qza \
	#--o-placements $OUT/insertion-placements.qza
## de novo
#qiime alignment mafft \
	#--i-sequences $OUT/classified_sequences/classification.qza \
	#--o-alignment $OUT/aligned-repset.qza
#qiime alignment mask \
	#--i-alignment $OUT/aligned-repset.qza \
	#--o-masked-alignment $OUT/masked-aligned-repset.qza
#qiime phylogeny fasttree \
	#--i-alignment $OUT/masked-aligned-repset.qza \
	#--o-tree $OUT/unrooted-tree.qza
#qiime phylogeny midpoint-root \
	#--i-tree $OUT/unrooted-tree.qza \
	#--o-rooted-tree $OUT/rooted-tree.qza

#Notes:
#primer
#如果引物，接头等长，可以用--p-trim 剪切；如果引物，接头不等长，可以使用cutadapt trim-paired #去除引物。

## 导入数据的可视化quality of sequences #9 mins
#time qiime demux summarize \
	#--i-data $OUT/paired-end-demux.qza \
	#--o-visualization $OUT/qual_viz.qzv
	
## check the quality file
#time qiime tools view $OUT/qual_viz.qzv
#1.哪个样本的测序深度最低？b14
#2.序列长度中位数是多少？250nts
#3.125位的中位数质量得分是多少？37

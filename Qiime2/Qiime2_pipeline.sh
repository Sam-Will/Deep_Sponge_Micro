## QIIME2 PIPELINE FOR BISECT 16S STUDY ###
# Sam Williams
# 05/12/2023

# First step is to do quality control on the raw illumina reads
# fastQC on the RawData from illumina
./fastqc -o /Users/sw17073/Desktop/16S_metagenomics/01_FastQC /Users/sw17073/Desktop/16S_metagenomics/00_Rawdata/*/*

# Activate the qiime2 conda environment
# Qiime2 is a comprehensive package for microbiome analysis
# Download here: https://docs.qiime2.org/2023.9/install/
conda activate qiime2-2023.7

# Next import the sequence data into Qiime2 as an qiime artifact
# Data is in CasavaOneEightSingleLanePerSampleDirFmt format

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path 00_rawdata/Pooled_samples \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path 02_QiimeArtifact/SMP-demux-paired-end.qza 

## So you can also visualise your data at this point. Firstly turning the .qza to a .qzv
qiime demux summarize \
  --i-data 02_QiimeArtifact/SMP-demux-paired-end.qza  \
  --o-visualization 02_QiimeArtifact/SMP-demux-paired-end.qzv

# Trim adaptors from the raw data with cutadapt
qiime cutadapt trim-paired \
--p-front-f CTGTCTCTTATACACATCTGTGYCAGCMGCCGCGGTAA \
--p-front-r CTGTCTCTTATACACATCTTTTGARTTTMYTTAACYGCC \
--i-demultiplexed-sequences 02_QiimeArtifact/SMP-demux-paired-end.qza \
--o-trimmed-sequences 03_Trimmed_seqs/SMP_demux_trimmed-seqs.qza \
--verbose | tee $PWD/Cutadapt_log.txt

# View post adaptor trimming
qiime demux summarize \
--i-data 03_Qiimecutadapt/SMP_demux_trimmed-seqs.qza \
--o-visualization 03_Qiimecutadapt/SMP_demux_trimmed-seqs.qzv

# Use DADA2 on denoise the data
# Choose your truncation lengths and trim based on quality plots
# You should have a minimum of 30 bp overlap between reads
qiime dada2 denoise-paired \
--i-demultiplexed-seqs 03_Qiimecutadapt/SMP_demux_trimmed-seqs.qza \
--p-trim-left-f 20 \
--p-trim-left-r 20 \
--p-trunc-len-f 269 \
--p-trunc-len-r 172 \
--p-n-threads 8 \
--o-table 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_table.qza \
--o-representative-sequences 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_rep-seqs.qza \
--o-denoising-stats 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_denoising-stats.qza 

# Download the Greengenes2 2022.10 full-length sequences classifier from here https://docs.qiime2.org/2023.7/data-resources/
# Used the full length over the 515F-805R as we have a longer sequence

qiime feature-classifier classify-sklearn \
  --i-classifier 05_classify/gg_2022_10_backbone_full_length.nb.qza \
  --i-reads 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_rep-seqs.qza \
  --o-classification 05_classify/SMP-truncr-172-truncf-269-trim20-DADA2_taxonomy.qza

# You can view a quick taxa bar blot with Qiime2
# You metadata file must be formatted corrected - See the docs
qiime taxa barplot \
  --i-table 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_table.qza \
  --i-taxonomy 05_classify/SMP-truncr-172-truncf-269-trim20-DADA2_taxonomy.qza \
  --m-metadata-file BISECT_metadata.tsv \
  --o-visualization 06_taxabarplot/SMP-truncr-172-truncf-269-trim20-DADA2_taxa-bar-plots.qsv
  
# Create a phylogenetic tree for the final piece of the phyloseq object
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_rep-seqs.qza \
  --o-alignment 07_phylo/aligned-rep-seqs.qza \
  --o-masked-alignment 07_phylo/masked-aligned-rep-seqs.qza \
  --o-tree 07_phylo/unrooted-tree.qza \
  --o-rooted-tree 07_phylo/rooted-tree.qza

# Export Qiime data ready for further analysis in phyloseq
qiime tools export --input-path 05_classify/SMP-truncr-172-truncf-269-trim20-DADA2_taxonomy.qza --output-path exported-taxonomy
qiime tools export --input-path 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_table.qza --output-path exported-feature-table
qiime tools export --input-path 04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_rep-seqs.qza --output-path exported-rep-seqs
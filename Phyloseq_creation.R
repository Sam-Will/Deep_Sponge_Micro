# An R script to produce the phyloseq object from the Qiime2 data processing

# Install various packages required for analysis
# Not required if already installed
Sys.setlocale("LC_ALL", "en_US.UTF-8")
install.packages("BiocManager")
BiocManager::install("GenomeInfoDb")
BiocManager::install("Biostrings")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("joey711/phyloseq")
install.packages("tidyverse")
install.packages("devtools")
devtools::install_github("jbisanz/qiime2R")
BiocManager::install("phyloseq")
BiocManager::install("decontam")

# Load the libraries
library(qiime2R)
library(phyloseq)
library(decontam)

# Set the working directory
setwd('/Users/saedwi/Qiime2')

# Create the phyloseq object for further analysis
BISECT.physeq<-qza_to_phyloseq(
  features="04_DADA2/SMP-truncr-172-truncf-269-trim20-DADA2_table.qza",
  tree="07_phylo/rooted-tree.qza",
  taxonomy="05_classify/SMP-truncr-172-truncf-269-trim20-DADA2_taxonomy.qza",
  metadata = "BISECT_metadata.tsv"
)

# Removing blanks, unwanted samples and decontam samples
# First lets remove the blank, gel extract and B00078
# B00078 was a test extraction and not used for the analysis
sponge_samples <- subset_samples(BISECT.physeq, Sample_or_Control == "True_sample")
sponge_samples <- subset_samples(BISECT.physeq, depth !=971) # Remove B00078
BISECT.physeq_pruned <- prune_samples(sample_names(sponge_samples), BISECT.physeq)

# Now we'll remove the 26 likely contaminant ASVs identified by decontam
# features taken from Decontam_combined_contam_frequency.csv in the 10_decontam folder
decon_features<- c("50633c4b837aa1a3b21742c0765270a4","ec1735efd303dce8c3a0043d92ff93c8","15157ac5595c109ea66f313a3eb5eb7a","f21dfad54cde3762dea1d00a41729507", "d838a3400b85af80da105b0545a6c789","d0bffe0ce5e1b758547fe2b490066b94","f1be6b0c6e4b3136ead03bc4d5c67c79","b9b990ccd54fd7edbe52563acac38b90","a04b27b43cec356541b84c9346d36125","63eed8bec991106ef19f974d5f2063f2", "c055aacdfae9d0169fc1c6f884292c4c","3464a7e5d5c33c0a88af22e879eb89a4","d2bad679382dfdf5b82f911661bf5dfb","dd606ddbf01aa42616b1135ac87f3fea","0d53eae08491290ff50a3c98788baf5b","5921c63ef5ace7407df7613b0113ea46","927ef3265e9e3355963d4408200e31f8")

# Get the list of all features and Identify the features you want to keep
all_features <- taxa_names(BISECT.physeq)
features_to_keep <- setdiff(all_features, decon_features)

# Prune the unwanted features from the phyloseq object
BISECT.physeq_decontam_pruned <- prune_taxa(features_to_keep, BISECT.physeq_pruned)

# Finally I want to update the sample names aka sponges (rows in phyloseq object)
# Currently they are not particularly meaningful so change to Sample_label which names the sponges by taxonomy
sample_metadata <- sample_data(BISECT.physeq_decontam_pruned)
sample_names(BISECT.physeq_decontam_pruned) <- sample_metadata$Sample_label
print(sample_names(BISECT.physeq_decontam_pruned))

# Tidy up
BISECT.physeq <- BISECT.physeq_decontam_pruned
rm(BISECT.physeq_decontam_pruned, BISECT.physeq_pruned, sponge_samples, decon_features, 
   all_features, features_to_keep, sample_metadata)

# Save the pruned phyloseq object
save(BISECT.physeq, file = "BISECT.physeq.RData")

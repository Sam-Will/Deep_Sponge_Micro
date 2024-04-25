# Deep_Sponge_Micro

Welcome to the **Deep_Sponge_Micro** GitHub repository. This repository is dedicated to hosting a collection of analysis scripts used for the research paper titled "Diversity and structure of the deep-sea sponge microbiome in the equatorial Atlantic Ocean."

## Repository Contents

- **Qiime2_pipeline.sh**: This script outlines the QIIME2 pipeline used for initial amplicon processing from raw reads into ASVs (Amplicon Sequence Variants) with assigned taxonomy.
- **Phyloseq_creation.R**: This R script is for creating the Phyloseq object from the Qiime2 pipeline output.
- **BISECT.physeq.RData**: This Rdata file contains the Phyloseq object used for analysis. You can load this into R and analyze the microbiome data with the phyloseq package.

## Data Availability

The sequence reads relevant to this research are deposited under the BioProject number: [PRJNA702029](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA702029).

---

*For more information or if you wish to contribute to this project, please feel free to open an issue or submit a pull request.*

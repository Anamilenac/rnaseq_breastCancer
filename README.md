# Characterize Gene Expression Differences Between Human Breast Cancer Subtypes

This repository contains code used for bioinformatics data analysis of the project " Characterize Gene Expression Differences Between Human Breast Cancer Subtypes"   
**Date:**  January 2024    
**Author:**  Ana Milena Castro Marquez  

# Data Availability

mRNA sequencing data used in the project is a subset of the data from Eswaran et al. 2012 (https://doi.org/10.1038/srep00264), obtained from the Gene Expression Omnibus (GEO) database (https://www.ncbi.nlm.nih.gov/geo) under accession GSE52194.

The dataset includes Illumina sequencing reads obtained through paired-end sequencing. Each sample is represented by two files, corresponding to read 1 and read 2.

The subset used in this project corresponds to three replicates for each of the following experimental groups, each originating from three different subtypes of human breast tumors, and three replicates of healthy control samples:
- **TNBC:** Triple-negative breast cancer
- **NonTNBC:** Non-triple-negative breast cancer
- **HER2:** HER2-positive breast cancer
- **Normal:** Healthy controls  

Replicates are identified with a number appended to the experimental group (e.g., HER21 for the first replicate in HER2), followed by reads _R1 and _R2.

# Data Analysis

The code is divided into 10 script files.

- 01_reads 
- 02_quality_control
- 05_indexing
- 06_1_list
- 06_mapping
- 07_samtobam
- 08_index_bam
- 09_feature_counts
- 10_clean_feature_counts
- 11_DEG

In order to replicate the code, names of analyses and output folders can be set.

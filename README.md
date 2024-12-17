# ProMix
Improving design and normalization of multiplex proteomics study

## 1. Introduction
The R script "promix_func.R" implements an iteration algorithm for improving design and normalization of multiplex proteomics study, described in Fang et al. 2024 (https://www.biorxiv.org/content/10.1101/2024.12.05.627093v1). It requires R and R package lme4. The following example has been tested under macOS Sonoma 14.4 with R version 4.4.2 and R packages lme4 v1.1-35.5 and Matrix v1.7-0. Please contact huatang@stanford.edu or hyfang@cnu.edu.cn for questions or bug report.


## 2. Usage
### 2.1 Prepare the input file with the experiment design and the raw protein abundance data

The input file includes the experiment design and the protein abundance data. The first three columns should be the sample (Samp), the subject (Subj), the run id (Run) and the raw protein abundance are followed from the fourth columns . The first line of the input file is the column names (Samp Subj Run Protein1_NAME Protein2_NAME ... ...) that they should be seperated by a space. One example file "input_promix.txt" is provided.

### 2.2 Run ProMix
Use "Rscript promix.R INPUT_FILE OUTPUT_FILE" in terminal for running ProMix. The output file will be in OUTPUT_FILE. For the example input file, the command in terminal can be "Rscript promix input_promix.txt output_promix.txt".

### 2.3 Explanation of output file
The output file includes the subject (the first column named Subj) and the protein abundance after normalization.


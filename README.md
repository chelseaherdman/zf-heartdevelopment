# zf-heartdevelopment

This repository assembles data and code from Herdman et al. 2020 zebrafish heart development time course.

## Contents

In this repository, you will find all analysis of the ribozero dataset from abundance estimates with kallisto to differential expression analysis with DESeq2.
This work is divided into three subfolders:

* kallisto_quant_combined_z11
    * Markdown file outlining our use of kallisto
    * Subfolders containing the log file, json and .tsv counts files for each sample
* Tximport
    * R script and Markdown file describing our use of tximport package
    * Bias corrected, gene summarized counts table produced by tximport
* DESeq2
    * R script and Markdown file
    * Diagnostic Figures subfolder
    * Figures subfolder
    * Differential expression analysis results table

## Experimental Materials & Methods

## Bioinformatic Analysis

### Authors

* **Chelsea Herdman** [Github](https://github.com/chelseaherdman)
* **Bradley L. Demarest** [Github](https://github.com/bdemarest)

Please contact Chelsea Herdman [cherdman@genetics.utah.edu] with any questions or concerns. 

***
### References
* Bray N, Pimentel H, Melsted P, Pachter L (2016), _Near-optimal probabilistic RNA-seq quantification_, Nature Biotechnology, 34, 525–527. doi:10.1038/nbt.3519  
* Soneson C, Love MI, Robinson MD (2015), _Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences_, F1000Research, 4, 1521. doi:10.12688/f1000research.7563.2  
* Love MI, Huber W, Anders S (2014), _Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2_, Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8

***


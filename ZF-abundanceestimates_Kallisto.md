
# Estimate Transcript Abundances using Kallisto

Add text here.

### Prepare kallisto index

Download fasta files for coding and non-coding annotated sequences from [ensembl](http://www.ensembl.org/info/data/ftp/index.html).

```
#! /bin/bash

kallisto index \
    -i GRCz11.r99.cdna.all.ncrna.kallisto_index \
    Danio_rerio.GRCz11.cdna.all.fa.gz Danio_rerio.GRCz11.ncrna.fa.gz
    
```

### Run kallisto in R using the package paralell

Run kallisto version 0.46.1 for ribozero total rnaseq timecourse experiment 
using GRCz11 and ensembl annotation release 99

```
library(data.table)
library(parallel)

# Load sample info sheet
sample_info = fread("../sample_info_14893R_20180221.txt")

sample_info = sample_info[, list(fastq_pairs=paste(paste(file.path("../fastq",
                                                                   read_1_filename),
                                                         file.path("../fastq",
                                                                   read_2_filename),
                                                   collapse=" "))), by=gnomex_id]

sample_info[, logfile_name:=paste("kallisto_quant", gnomex_id, "log.txt", sep="_")]

command_vec = paste(
    "mkdir -p", sample_info$gnomex_id, ";",
    "kallisto quant",
    "-i /data3/cherdman/genome_z11/GRCz11.r99.cdna.all.ncrna.kallisto_index",
    "-o", sample_info$gnomex_id,
    "-b 100",
    "-t 6",
    sample_info$fastq_pairs,
    ">>", file.path(".", sample_info$gnomex_id, sample_info$logfile_name), "2>&1")

cl = makeCluster(4)

```
### Made long-form datatable of kallisto results in R

```
library(data.table)

file_name_vec = list.files(path=".",
                           pattern="abundance.tsv",
                           recursive=TRUE,
                           full.names=TRUE)

gnomex_id_vec = basename(dirname(file_name_vec))

results_list = list()

for (j in seq_along(file_name_vec)) {
  file_name = file_name_vec[j]
  gnomex_id = gnomex_id_vec[j]

  tmp = fread(file_name)
  set(tmp, j="gnomex_id", value=gnomex_id)
  results_list[[j]] = tmp
}

ktab = rbindlist(results_list)
setnames(ktab, old="target_id", new="ensembl_transcript_id")
ktab[, ensembl_transcript_id:=sapply(
          strsplit(ensembl_transcript_id, 
          split=".", fixed=TRUE), 
    function(x) x[1])]

fwrite(ktab, file="kallisto_quant_z11_14893R_n22_long_form.txt", sep="\t")

```

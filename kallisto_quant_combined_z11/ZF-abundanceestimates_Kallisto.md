
# Estimate Transcript Abundances using Kallisto

### Prepare kallisto index

Download fasta files for coding and non-coding annotated sequences from [ensembl](http://www.ensembl.org/info/data/ftp/index.html) to a directory named
"genome".

Run the index function in a bash command line from within the "genome" directory.

```
kallisto index \
    -i GRCz11.r99.cdna.all.ncrna.kallisto_index \
    Danio_rerio.GRCz11.cdna.all.fa.gz Danio_rerio.GRCz11.ncrna.fa.gz
```
### Download fastq files from this dataset

Download fastq files from [CvDC datahub](https://b2b.hci.utah.edu/gnomex/), experiment id 14893R and store in a folder named "fastq" in the same parent directory as "genome".

### Run kallisto in R using the package paralell

Run kallisto version 0.46.1 for ribozero total rnaseq timecourse experiment 
using GRCz11 and ensembl annotation release 99 in a directory named 
kallisto_quant_combined_z11.

```
mkdir kallisto_quant_combined_z11
cd kallisto_quant_combined_z11
```
Start R.

```
library(data.table)
library(parallel)

# Load sample info sheet and collapse fastq paired filenames by sample id (gnomex_id in our case).

sample_info = fread("kallisto_sample_info_14893R.txt")

sample_info = sample_info[, list(fastq_pairs=paste(paste(file.path("../fastq",
                                                                   read_1_filename),
                                                         file.path("../fastq",
                                                                   read_2_filename),
                                                   collapse=" "))), by=gnomex_id]

sample_info[, logfile_name:=paste("kallisto_quant", gnomex_id, "log.txt", sep="_")]

command_vec = paste(
    "mkdir -p", sample_info$gnomex_id, ";",
    "kallisto quant",
    "-i ../genome/GRCz11.r99.cdna.all.ncrna.kallisto_index",
    "-o", sample_info$gnomex_id,
    "-b 100",
    "-t 6",
    sample_info$fastq_pairs,
    ">>", file.path(".", sample_info$gnomex_id, sample_info$logfile_name), "2>&1")

cl = makeCluster(4)

```

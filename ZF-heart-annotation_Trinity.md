
# Zebrafish Heart Annotation using Trinity

We were interested in verifying heart specific annotations of transcripts, and especially 3'UTRs using the high resolution transcriptome sequenced in this project.

## Getting Started

In order to build a heart specific annotation of total RNA transcripts during development, we employed Trinity. 

### Installation

Download Trinity.

```
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.9.1/trinityrnaseq-v2.9.1.FULL.tar.gz
```
Trinity requires cmake3

```
yum install cmake3
```
Trinity also requires:
* bowtie2
* jellyfish
* salmon
* python2.7 or 3 with numpy

Compile Trinity and run test by following instructions here:

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity

## Ran Trinity on our samples

### Prepped sample info table for Trinity in R

Trinity will take a tab-delimited text file with four columns:
* condition
* sample/replicate id
* read 1 file name
* read 2 file name

```
library(data.table)

tab = fread("../sample_info_14893R_20180221.txt")

tab2 = data.table(time_point=tab$time_point,
                  rep_id=gsub("^(\\d)_.+$", "\\1",  tab$sample_id),
                  lane_id=tab$lane, rep_lane="",
                  read1=file.path("/data3/cherdman/heart_timecourse/ribozero_rnaseq/fastq", 
                                  tab$read_1_filename),
                  read2=file.path("/data3/cherdman/heart_timecourse /ribozero_rnaseq/fastq", 
                                  tab$read_2_filename))
```
Make a unique fastq pair identifier by combining lane and replicated id.
Be careful that Trinity doesn't treat these all as biological replicates.

```
tab2[, rep_lane:=paste(time_point, "_rep", rep_id, "_lane", lane_id, sep="")]

set(tab2, j="rep_id", value=NULL)
set(tab2, j="lane_id", value=NULL)

fwrite(tab2, file="trinity_sample_info_20200115.txt", sep="\t", col.names=FALSE)

```

### Ran Trinity in shell script

Ran Trinity v2.9.0 using chosen memory and ram requirements specific to our server.
Printed both standard out and standard error to log file. 

```
Trinity \
  --seqType fq \
  --max_memory 256G \
  --CPU 32 \
  --samples_file trinity_sample_info_20200115.txt \
  --SS_lib_type FR > trinity_log_20200123.txt 2>&1
```

## Authors

* **Chelsea Herdman** [Github](https://github.com/chelseaherdman)
* **Bradley L. Demarest** [Github](https://github.com/bdemarest)

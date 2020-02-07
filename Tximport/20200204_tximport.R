#' ---
#' title: "Import transcript abundances from kallisto using tximport"
#' author: "Chelsea Herdman"
#' date: "February 5th, 2020"
#' output: github_document
#' ---
#'
#' In order to perform differential expression analysis using DESeq2 on the 
#' estimated transcript abundances produced by kallisto, we used the tximport
#' package in R.
#'
#' We followed the vignette provided by Michael I. Love, Charlotte Soneson and 
#' Mark D. Robinson found [here](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon_with_inferential_replicates)
#' 
#' #### Prepare transcript to gene dataframe for gene-level summarization
#' 
#' Load required libraries.
#' 
#+ libraries, message=FALSE, error=FALSE, warning=FALSE
library(tximport)
library(data.table)
library(ensembldb)
library(RMariaDB)
library(here)

#' Fetch gene annotations from Ensembl and create one to one 
#' transcript id/gene id data frame.

txdb = makeTxDbFromEnsembl(organism="Danio rerio",
                           release=99,
                           circ_seqs=DEFAULT_CIRC_SEQS,
                           server="ensembldb.ensembl.org",
                           username="anonymous", password=NULL, port=0L,
                           tx_attrib=NULL)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

length(unique(tx2gene$TXNAME))

length(unique(tx2gene$GENEID))

#' #### Prepare named vector of files
#' Use kallisto abundance.tsv files. These are listed in alphabetical order in the 
#' kallisto_out directory so ensure naming associates accurately and that names 
#' associate with sample info table that will be used for DESeq2.

file_path_vec = list.files(path=here("kallisto_quant_combined_z11"),
                   pattern="abundance.tsv",
                   recursive = TRUE,
                   full.names = TRUE)
all(file.exists(file_path_vec))

file_tab = data.table(file_path=file_path_vec,
                      gnomex_id=basename(dirname(file_path_vec)))
file_tab[, sort_order:=as.integer(gsub("^14893X", "", gnomex_id))]
setorder(file_tab, sort_order)

files = file_tab$file_path
names(files) = file_tab$gnomex_id

#'
#' #### Run tximport for gene level estimation
#'
#' DESeq2 will be run both using the *original counts and offset* method and the 
#' *bias corrected counts without an offset* method in order to compare the two 
#' possibilities of correcting for transcript length. 
#'
#' ##### Original counts and offset method
#' This txi object can be used directly with downstream DESeq2 analyses using the
#' DESeqDataSetFromTximport function.
#' 
txi.kallisto.offset = tximport(files, type= "kallisto", tx2gene = tx2gene, 
                               ignoreTxVersion = TRUE)
txi.kallisto.offset$counts[1:6, 1:6]
txi.kallisto.offset$abundance[1:6, 1:6]

#' ##### Bias corrected counts without an offset
txi.kallisto.biascorr = tximport(files, type= "kallisto", tx2gene = tx2gene, 
                                 ignoreTxVersion = TRUE, 
                                 countsFromAbundance="lengthScaledTPM")
txi.kallisto.biascorr$counts[1:6, 1:6]

counts_tab = as.data.table(txi.kallisto.biascorr$counts)
counts_tab$ensembl_gene_id = rownames(txi.kallisto.biascorr$counts)
setcolorder(counts_tab, c("ensembl_gene_id", paste("14893X", 1:22, sep="")))

#' Save the counts table in order to use as input for DESeq2 (*see DESeq2 folder*)
fwrite(counts_tab, file=here("Tximport", "20200204_ribozero_counts_fromtximport_biascorrected.txt.gz"), sep="\t")


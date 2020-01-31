Tx Import for DESeq2
================
Chelsea Herdman
30/01/2020

In order to perform differential expression analysis using DESeq2 on the
estimated transcript abundances produced by kallisto, we used the
tximport package in R.

We followed the vignette provided by Michael I. Love, Charlotte Soneson
and Mark D. Robinson found
[here](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon_with_inferential_replicates)

## Prepare transcript to gene dataframe for gene-level summarization

Load required libraries.

``` r
library(tximport)
library(data.table)
library(ensembldb)
library(RMariaDB)
```

Fetch gene annotations from Ensembl and create one to one transcript
id/gene id data frame.

``` r
txdb = makeTxDbFromEnsembl(organism="Danio rerio",
                    release=99,
                    circ_seqs=DEFAULT_CIRC_SEQS,
                    server="ensembldb.ensembl.org",
                    username="anonymous", password=NULL, port=0L,
                    tx_attrib=NULL)
```

    ## Fetch transcripts and genes from Ensembl ... OK
    ## Fetch exons and CDS from Ensembl ... OK
    ## Fetch chromosome names and lengths from Ensembl ...OK
    ## Gather the metadata ... OK
    ## Make the TxDb object ... OK

``` r
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
head(tx2gene)
```

    ##               TXNAME             GENEID
    ## 1 ENSDART00000193350 ENSDARG00000114632
    ## 2 ENSDART00000181299 ENSDARG00000109882
    ## 3 ENSDART00000180768 ENSDARG00000111634
    ## 4 ENSDART00000184307 ENSDARG00000114112
    ## 5 ENSDART00000186443 ENSDARG00000116981
    ## 6 ENSDART00000189953 ENSDARG00000109637

``` r
length(unique(tx2gene$TXNAME))
```

    ## [1] 65905

``` r
length(unique(tx2gene$GENEID))
```

    ## [1] 37241

## Prepare named vector of files

Use kallisto abundance.h5 files. These are listed in alphabetical order
in the kallisto\_out directory so ensure naming associates accurately
and that names associate with sample info table that will be used for
DESeq2.

``` r
files = list.files(path="../../kallisto_out/",
                   pattern="abundance.h5",
                   recursive = TRUE,
                   full.names = TRUE)

all(file.exists(files))
```

    ## [1] TRUE

``` r
#[1] TRUE

names(files) = paste0("14893X", c(1, 10:19, 2, 20:22, 3:9))

files
```

    ##                                     14893X1 
    ##  "../../kallisto_out//14893X1/abundance.h5" 
    ##                                    14893X10 
    ## "../../kallisto_out//14893X10/abundance.h5" 
    ##                                    14893X11 
    ## "../../kallisto_out//14893X11/abundance.h5" 
    ##                                    14893X12 
    ## "../../kallisto_out//14893X12/abundance.h5" 
    ##                                    14893X13 
    ## "../../kallisto_out//14893X13/abundance.h5" 
    ##                                    14893X14 
    ## "../../kallisto_out//14893X14/abundance.h5" 
    ##                                    14893X15 
    ## "../../kallisto_out//14893X15/abundance.h5" 
    ##                                    14893X16 
    ## "../../kallisto_out//14893X16/abundance.h5" 
    ##                                    14893X17 
    ## "../../kallisto_out//14893X17/abundance.h5" 
    ##                                    14893X18 
    ## "../../kallisto_out//14893X18/abundance.h5" 
    ##                                    14893X19 
    ## "../../kallisto_out//14893X19/abundance.h5" 
    ##                                     14893X2 
    ##  "../../kallisto_out//14893X2/abundance.h5" 
    ##                                    14893X20 
    ## "../../kallisto_out//14893X20/abundance.h5" 
    ##                                    14893X21 
    ## "../../kallisto_out//14893X21/abundance.h5" 
    ##                                    14893X22 
    ## "../../kallisto_out//14893X22/abundance.h5" 
    ##                                     14893X3 
    ##  "../../kallisto_out//14893X3/abundance.h5" 
    ##                                     14893X4 
    ##  "../../kallisto_out//14893X4/abundance.h5" 
    ##                                     14893X5 
    ##  "../../kallisto_out//14893X5/abundance.h5" 
    ##                                     14893X6 
    ##  "../../kallisto_out//14893X6/abundance.h5" 
    ##                                     14893X7 
    ##  "../../kallisto_out//14893X7/abundance.h5" 
    ##                                     14893X8 
    ##  "../../kallisto_out//14893X8/abundance.h5" 
    ##                                     14893X9 
    ##  "../../kallisto_out//14893X9/abundance.h5"

## Run tximport for gene level estimation

``` r
txi.kallisto = tximport(files, type= "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
```

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length
    ## summarizing inferential replicates

-----

## References

  - Bray N, Pimentel H, Melsted P, Pachter L (2016), *Near-optimal
    probabilistic RNA-seq quantification*, Nature Biotechnology, 34,
    525–527. <doi:10.1038/nbt.3519>  
  - Soneson C, Love MI, Robinson MD (2015), *Differential analyses for
    RNA-seq: transcript-level estimates improve gene-level inferences*,
    F1000Research, 4, 1521. <doi:10.12688/f1000research.7563.2>  
  - Love MI, Huber W, Anders S (2014), *Moderated estimation of fold
    change and dispersion for RNA-seq data with DESeq2*, Genome Biology,
    15, 550. <doi:10.1186/s13059-014-0550-8>

-----

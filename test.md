DESeq2 on RiboZero (total RNA) Data
================

## Load kallisto results table and annotate

Load required packages, read in results table and sample info table.
Fetch gene annotations from Biomart using GRCz11, annotation version 99.

``` r
library(data.table)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(ggrepel)
```

``` r
ktab = fread("../../kallisto_out/kallisto_quant_z11_14893R_n22_long_form.txt")

head(ktab)
```

    ##    ensembl_transcript_id length eff_length est_counts      tpm gnomex_id
    ## 1:    ENSDART00000189226     10    6.33333          0 0.000000   14893X1
    ## 2:    ENSDART00000189431     11    5.75000          0 0.000000   14893X1
    ## 3:    ENSDART00000160762    345  156.98900          2 0.123762   14893X1
    ## 4:    ENSDART00000170804    366  176.21400          8 0.441036   14893X1
    ## 5:    ENSDART00000165410    350  161.50100          0 0.000000   14893X1
    ## 6:    ENSDART00000163675    339  151.55900          0 0.000000   14893X1

``` r
sample_info = fread("../../ribozero_sample_info.txt")

head(sample_info)
```

    ##    gnomex_id rep_time time_point
    ## 1:   14893X1  1_24hpf      24hpf
    ## 2:   14893X2  1_36hpf      36hpf
    ## 3:   14893X3  1_48hpf      48hpf
    ## 4:   14893X4  1_60hpf      60hpf
    ## 5:   14893X5  1_72hpf      72hpf
    ## 6:   14893X6  2_24hpf      24hpf

``` r
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
               dataset="drerio_gene_ensembl",
               host="jan2020.archive.ensembl.org")

annot = getBM(mart=mart,
              attributes=c("ensembl_gene_id",
                           "ensembl_transcript_id",
                           "external_gene_name",
                           "chromosome_name",
                           "gene_biotype"))

annot = as.data.table(annot)
```

Merge in gene annotation information and sum estimated counts on a
gene-by-gene
basis.

``` r
ktab = merge(ktab, annot[, list(ensembl_gene_id, ensembl_transcript_id)])

dim(ktab)
```

    ## [1] 1449580       7

``` r
length(unique(ktab$ensembl_transcript_id))
```

    ## [1] 65890

``` r
length(unique(ktab$ensembl_gene_id))
```

    ## [1] 37227

``` r
wide_tab = dcast(ktab, ensembl_gene_id ~ gnomex_id, value.var="est_counts",
                 fun.aggregate=sum)

setcolorder(wide_tab, c("ensembl_gene_id", sample_info$gnomex_id))
```

DESeq2 requires a counts matrix of integers.

``` r
counts = as.matrix(wide_tab[, -"ensembl_gene_id"])
rownames(counts) = wide_tab$ensembl_gene_id

storage.mode(counts) = "integer"

counts[1:6, 1:6]
```

    ##                    14893X1 14893X2 14893X3 14893X4 14893X5 14893X6
    ## ENSDARG00000000001     374     591     654     515     659     382
    ## ENSDARG00000000002    1199    1739    2833    2851    2885     958
    ## ENSDARG00000000018    4476    5743    4855    2571    3121    2772
    ## ENSDARG00000000019    3346    2919    3523    3638    3963    2973
    ## ENSDARG00000000068    6639     904     734    1118    1898    5708
    ## ENSDARG00000000069    1313    1175    1344    1233    1490     856

``` r
str(counts)
```

    ##  int [1:37227, 1:22] 374 1199 4476 3346 6639 1313 4539 11915 590 365 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:37227] "ENSDARG00000000001" "ENSDARG00000000002" "ENSDARG00000000018" "ENSDARG00000000019" ...
    ##   ..$ : chr [1:22] "14893X1" "14893X2" "14893X3" "14893X4" ...

## Run DESeq2

Load required packages, read in results table and sample info table.
Fetch gene annotations from Biomart using GRCz11, annotation version 99.

``` r
ggplot(mtcars, aes(x=mpg, y=hp)) + geom_point()
```

![](test_files/figure-gfm/prep%20cData-1.png)<!-- -->

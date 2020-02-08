# Create a tab-delimited text annotation file,
# to prevent the delays caused by repeatedly running getBM(), etc.

library(biomaRt)
library(data.table)

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
              dataset="drerio_gene_ensembl",
              host="jan2020.archive.ensembl.org")

annot = as.data.table(getBM(mart=mart,
                             attributes=c("ensembl_gene_id",
                                          "external_gene_name",
                                          "chromosome_name",
                                          "gene_biotype")))


fwrite(annot, file=here("DESeq2", "gene_annot_biomart_z11_r99.txt.gz"), sep="\t")

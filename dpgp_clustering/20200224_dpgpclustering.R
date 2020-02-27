#' ---
#' title: "Import clustering results from DPGP output"
#' author: "Chelsea Herdman"
#' date: "February 24th, 2020"
#' output: github_document
#' ---
#'
#' In order to perform GO term analysis on the clustering results from DPGP
#' we imported optimal clustering information from the DPGP run and created
#' gene lists per cluster
#'
#' #### Evaluate number of genes per cluster and probability statistic
#' 
#' Load required libraries.
#' 
#+ libraries, message=FALSE, error=FALSE, warning=FALSE
library(data.table)
library(here)
library(ggplot2)

#' Import optimal clustering table

dpgp_res = fread(here("dpgp_clustering", "dpgp_default_vsd_fdr5pct_replicates_20200222_optimal_clustering.txt"))

#' Make summary table of cluster information
sum_res = dpgp_res[, list(n_genes=length(gene),
                          mean_prob=mean(probability)),
                   by="cluster"]

#' Plot probability results for each cluster

#+ density plots, width=6, height=4

p1 = ggplot(dpgp_res, aes(x=probability)) +
        theme_bw() +
        geom_density(fill="#7bccc4") +
        scale_x_continuous(name= "Probability") +
        scale_y_continuous(name= "Density") +
        ggtitle("Density plot of cluster assignment probability") +
        facet_wrap(~ cluster)
       
cluster_vec = unique(dpgp_res$cluster) 

res_list = list()

for (i in seq_along(cluster_vec)) {
        tmp = dpgp_res[cluster == i]
        tmp = tmp[, list(gene)]
        fwrite(tmp, file=here("dpgp_clustering", "cluster_genelists", 
                              paste("cluster", i, ".txt", sep="")), 
               col.names = FALSE)
}

#' Use the cluster gene lists as input on pantherdb.org and the background
#' list of 35562 genes (subsetted gene list that was used as input for DESeq2)
#' 
#' 

#===============================================================================
# Read in panther results tables.

file_tab = data.table(filename=list.files(path=here("dpgp_clustering", "panther_output"), 
                                          pattern="FDR_res.txt", 
                                          full.names=TRUE))

file_tab[, cluster_id:=gsub(pattern="(cluster\\d{1,2})_.+$", 
                            replacement="\\1",
                            x=basename(filename))]

# Example of original column names.
# [1] "GO biological process complete"  "backgroundlist_35562genes.txt - REFLIST (23232)"  
# [3] "cluster1.txt (849)"             "cluster1.txt (expected)"       
# [5] "cluster1.txt (over/under)"      "cluster1.txt (fold Enrichment)"
# [7] "cluster1.txt (raw P-value)"     "cluster1.txt (FDR)" 


new_col_names = c("go_term", "reflist_count", "cluster_count",
                  "cluster_expected", "enrichment_direction", "fold_enrichment",
                  "pvalue", "padj")

res_list = list()

for (i in seq(nrow(file_tab))) {
        tmp = fread(file_tab[i, filename], skip=11)
        
        cluster_size = gsub(pattern="cluster\\d{1,2}\\.txt \\((\\d{1,5})\\)$", 
                            replacement="\\1", 
                            x=names(tmp)[3])
        
        setnames(tmp, new_col_names)
        
        set(tmp, j="cluster_id", value=file_tab[i, cluster_id])
        set(tmp, j="cluster_size", value=cluster_size)
        
        res_list[[i]] = tmp
}

tab = rbindlist(res_list)

tab[, go_term_id:=gsub(pattern=".+\\((GO:\\d{1,7})\\)$", 
                       replacement="\\1",
                       x=go_term)]

#-------------------------------------------------------------------------------
# What is going on with enrichment results set to "< 0.01"?

tail(sort(table(tab$fold_enrichment)))
table(tab[fold_enrichment == "< 0.01", enrichment_direction])

# GO terms of rows with "< 0.01" all have "-" enrichment_direction.

xtab = tab[fold_enrichment == "< 0.01" & enrichment_direction == "-"]
table(xtab$cluster_count)
head((table(xtab$reflist_count)))

# Cluster counts are all zeros, and reflist counts are very small.

# Check if it makes sense to remove GO terms with "-" and "< 0.01".
xtab2 = tab[go_term_id %in% unique(xtab$go_term_id)]
sum(xtab$padj < 0.05)
sum(xtab2$padj < 0.05)
length(unique(xtab$go_term_id))
length(unique(tab$go_term_id))

# Conclusion:
# No, we should not remove these because most of the affected GO terms
# have FDR < 0.01 in one or more clusters.

#-------------------------------------------------------------------------------
# Fix fold_enrichment column.
# It is class 'character' because it has (very) many "< 0.01" values.
# Here, we convert those to numeric 0.01.
tab[, fold_enrichment:=as.numeric(gsub(pattern="< ", 
                                       replacement="", 
                                       fold_enrichment))]


#fwrite(tab, here("dpgp_clustering", "panther_output", "20200226_combinedpantherresultscluster1-7.txt"), sep="\t")

# Diagnostic plots for panther results table.
h1 = ggplot(tab, aes(x=pvalue)) +
     geom_histogram(bins=200, color="grey30", fill="grey80")

h2 = ggplot(tab, aes(x=-log10(pvalue))) +
        geom_histogram(bins=200, color="grey30", fill="grey80")

length(unique(tab[padj < 0.01, go_term_id]))
# [1] 654

#-------------------------------------------------------------------------------
# Select for plotting: GO terms which have padj > 0.01 in at least one cluster.

subtab = tab[go_term_id %in% unique(tab[padj < 0.01, go_term_id])]



# To do:
# (1) Change fold_enrichment such that they are still log2 values, but bars
#     extend in same direction for both "-" and "+".
# 
# (3) hclust of rows (GO terms). Need matrix of fold_enrichment values
#     where columns are cluster_id and rows are go_term_id.
# (4) Crazy idea. Interactive plot where user hovers over bar to see
#     cluster_count, ref count, maybe gene list?
# (5) Bring GO term and ensembl_gene_ids from biomart,
#     to make the gene ids available for interpreting GO enrichment results.


# Set "< 0.01" enrichment values to NA
subtab[fold_enrichment <= 0.01, fold_enrichment:=NA_real_]

subtab[, log2_fold_enrichment:=log2(fold_enrichment)]
subtab[, abs_log2_fold_enrichment:=abs(log2_fold_enrichment)]
subtab[, fdr_group:=ifelse(padj < 0.01, "signif", "nonsignif")]

# Use hclust to set y-axis order.
fe_tab = dcast(data=subtab, 
               formula=go_term_id ~ cluster_id, 
               value.var="log2_fold_enrichment")

fe_mat = as.matrix(fe_tab, rownames="go_term_id")
fe_mat[is.na(fe_mat)] <- 0.0
d = dist(fe_mat, method="euclidean")
cl = hclust(d, method="average")
pdf(here("dpgp_clustering", "hlust_test_euc_average.pdf"), width=32, height=16)
plot(cl)
dev.off()

go_tab = data.table(go_term_id=cl$labels[cl$order])
go_tab$go_term_order = seq(nrow(go_tab))

go_tab = merge(x=go_tab, 
               y=subtab[, list(go_term_id, go_term)][!duplicated(subtab[, list(go_term_id, go_term)])], 
               all.x=TRUE, 
               all.y=FALSE)
setorder(go_tab, go_term_order)

subtab[, go_term:=factor(go_term, levels=go_tab$go_term)]

# Set four-level factor for enrichment direction X fdr group.

subtab$fdr_enrichment_group = NA_character_

subtab[enrichment_direction == "+" & fdr_group == "signif", 
       fdr_enrichment_group:="Enriched - FDR < 0.01"]

subtab[enrichment_direction == "+" & fdr_group == "nonsignif", 
       fdr_enrichment_group:="Enriched - FDR >= 0.01"]

subtab[enrichment_direction == "-" & fdr_group == "signif", 
       fdr_enrichment_group:="Depleted - FDR < 0.01"]

subtab[enrichment_direction == "-" & fdr_group == "nonsignif", 
       fdr_enrichment_group:="Depleted - FDR >= 0.01"]


fdr_enrichment_group_colors_2 = c(
        "Depleted - FDR < 0.01"="#e08214", ## dark orange
        "Depleted - FDR >= 0.01"="#fee0b6", ## lighter orange
        "Enriched - FDR >= 0.01"="#d8daeb", ## lighter purple
        "Enriched - FDR < 0.01"="#8073ac" ## dark purple
)



s1 = ggplot(subtab, aes(y=abs_log2_fold_enrichment, x=go_term, 
                        fill=fdr_enrichment_group)) +
        theme_bw() +
        geom_bar(colour="grey30", stat="identity", 
                 size=0.2, show.legend=FALSE) +
        scale_fill_manual(values=fdr_enrichment_group_colors_2) +
        coord_flip() +
        facet_grid(. ~ cluster_id)

ggsave(here("dpgp_clustering", "20200226_long_test_GOterm_plots.pdf"), plot=s1, 
       height=72, width=14, limitsize=FALSE)


#-------------------------------------------------------------------------------
# Dig into how to interpret and plot the fold_enrichment "< 0.01" values.
sx_tab = tab[, list(reflist_count=reflist_count[1], 
                    reflist_check=all(reflist_count == reflist_count[1])), 
             by=go_term_id]

tail(sort(table(sx_tab$reflist_count)))









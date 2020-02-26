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

dpgp_res = fread(here("dpgp_clustering/dpgp_default_vsd_fdr5pct_replicates_20200222_optimal_clustering.txt"))

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

file_tab = data.table(filename=list.files(path=here("dpgp_clustering", "panther_output"), 
                                          pattern="FDR_res.txt", 
                                          full.names=TRUE))

file_tab[, cluster_id:=gsub(pattern="cluster(\\d{1,2})_.+$", 
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

#fwrite(tab, here("dpgp_clustering", "panther_output", "20200226_combinedpantherresultscluster1-7.txt"), sep="\t")




S1 = ggplot(revigo_res, aes(x=cluster_id, y=go_term, 
                            size=fold_enrichment, color=padj)) + 
        geom_point(alpha = 0.8) + 
        theme_classic() +
        scale_fill_gradient(low = "red2",  high = "mediumblue", space = "Lab")
#   scale_size(range = c(2, 8))

# ggsave("20200120_testGOtermsplots.png", plot = S1, height=12, width=12, units="in", dpi = 150)

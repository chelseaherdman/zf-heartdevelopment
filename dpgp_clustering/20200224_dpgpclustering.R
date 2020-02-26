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


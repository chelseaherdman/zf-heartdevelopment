




sampleTable = fread("DESeq2/ribozero_sample_info.txt")

cData = data.frame(time_point=factor(sampleTable$time_point,
                                     levels=c("24hpf", "36hpf", "48hpf", "60hpf",
                                              "72hpf")),
                   replicate_id=factor(sampleTable$rep_id,
                                       levels=paste("rep_", 1:5, sep="")))

rownames(cData) = sampleTable$gnomex_id

# dds <- DESeqDataSetFromTximport(txi.kallisto, cData, ~time_point)
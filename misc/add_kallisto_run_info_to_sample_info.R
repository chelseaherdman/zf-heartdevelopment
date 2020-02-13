# From kallisto quant report, get total reads and assigned reads using rjson.
# Add to sample info sheet.

library(data.table)
library(rjson)
library(here)


json_file_vec = list.files(path=here("kallisto_quant_combined_z11"),
                           pattern="run_info.json", full.names=TRUE, recursive=TRUE)

gnomex_id_vec = basename(dirname(json_file_vec))

res_list = list()

for (i in seq_along(json_file_vec)) {
  tmp_json = fromJSON(file=json_file_vec[i])
  tmp_tab  = as.data.table(tmp_json)
  tmp_tab[, gnomex_id:=gnomex_id_vec[i]]
  res_list[[i]] = tmp_tab
}

jtab = rbindlist(res_list)


sample_info = fread(here("DESeq2", "ribozero_sample_info.txt"))

new_sample_info = merge(sample_info,
                        jtab[, list(n_processed, n_pseudoaligned, n_unique,
                                    p_pseudoaligned, p_unique, gnomex_id)],
                        by="gnomex_id")

# Make sure the table is sorted by gnomex id X1-X22.
setorder(new_sample_info, rep_time)

fwrite(new_sample_info, file=here("DESeq2", "extended_sample_info_14893R.txt"), sep="\t")




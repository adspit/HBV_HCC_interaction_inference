# prepare file for network propagation

LIHC_mutations_ppr_file_hbv = results_df
dim(LIHC_mutations_ppr_file_hbv)
# [1] 8765   14
colnames(LIHC_mutations_ppr_file_hbv)[1] =  "TCGA_gene_name"

# merge TCGA differential mutation analysis and HBV interactions list
ppr_file_hbv = merge(selected_preys_hbv, LIHC_mutations_ppr_file_hbv, by = "TCGA_gene_name", all = T)
dim(ppr_file_hbv)
# 9216    19

write.table(ppr_file_hbv, "network_propagation/propagation_file_LIHC_TCGA_HBV.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

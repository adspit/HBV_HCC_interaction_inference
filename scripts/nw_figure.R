#########################################
###### network figure HBV extended cohort
#########################################

# network figure:
viral_human_interactions = MIST_file_HBV[,c("TCGA_gene_name","Bait", 'MIST')]
dim(viral_human_interactions)
# [1] 3831    3


significant_human_protein = dat10[dat10$combined_prop_fdr_2<0.1,]
colnames(significant_human_protein)[1]="TCGA_gene_name"
dim(significant_human_protein)
# [1] 53  21



# select the genes in the viral_human nw that were found to be significant after the network propation
viral_human_interactions = merge(viral_human_interactions[viral_human_interactions$TCGA_gene_name%in%significant_human_protein$TCGA_gene_name,],significant_human_protein[,c("TCGA_gene_name","point_col", "Deviance_hepB")],by="TCGA_gene_name")
dim(viral_human_interactions)
# [1] 163   5

viral_human_interactions$Deviance_hepB[viral_human_interactions$point_col%in%"deletion"] = - viral_human_interactions$Deviance_hepB[viral_human_interactions$point_col%in%"deletion"]

dim(viral_human_interactions)
# [1] 163 5
viral_human_interactions = viral_human_interactions[viral_human_interactions$point_col!='non_significant',]
dim(viral_human_interactions)
# [1] 163   5

#write.table(viral_human_interactions, "HBV_viral_human_interactions.txt", sep="\t", col.names = T, row.names = F, quote = F)
write.table(viral_human_interactions, "HBV_viral_human_interactions_fdr_10.txt", sep="\t", col.names = T, row.names = F, quote = F)




# high confidence interaction list network based on MIST/Compass scores
HIC_viral_human = MIST_file_HBV[MIST_file_HBV$MIST>0.75,c("TCGA_gene_name","Bait",'MIST')]
dim(HIC_viral_human)
# [1] 145   3



# high confidence interaction list network that contains significant genes
HIC_viral_human_interactions = merge(HIC_viral_human[HIC_viral_human$TCGA_gene_name%in%significant_human_protein$TCGA_gene_name,],significant_human_protein[,c("TCGA_gene_name","point_col", "Deviance_hepB"), ], by="TCGA_gene_name")
dim(HIC_viral_human_interactions)
# [1] 14 5

#HIC_viral_human_interactions = HIC_viral_human_interactions[HIC_viral_human_interactions$point_col!='non_significant',]
#HIC_viral_human_interactions$Deviance_hepB[HIC_viral_human_interactions$point_col%in%"decreased_HBV_mut"] = - HIC_viral_human_interactions$Deviance_hepB[HIC_viral_human_interactions$point_col%in%"decreased_HBV_mut"]
write.table(HIC_viral_human_interactions, "HIC_HBV_viral_human_interactions_fdr_10_MIST_75.txt", sep="\t", col.names = T, row.names = F, quote = F)


### until here, we have the networks from our data - without the Reactome human PPI nw




# read in the human Reactome PPI nw
generic_human_PPI = read.table("../updated_FIsInGene_031516_with_annotations.txt", h=T, sep="\t", stringsAsFactors = F)
dim(generic_human_PPI)
# [1] 229300      5


# select only genes that are found in the significant list after the propagation algorithm
reduced_generic_PPI = generic_human_PPI[generic_human_PPI$Gene1%in%significant_human_protein$TCGA_gene_name & generic_human_PPI$Gene2%in% significant_human_protein$TCGA_gene_name,1:2 ]
colnames(reduced_generic_PPI)[1] = "TCGA_gene_name"

dim(reduced_generic_PPI)
# [1] 89  2
# add column for increased/decreased mutation rate
all.equal.character(reduced_generic_PPI$TCGA_gene_name, as.character(significant_human_protein$TCGA_gene_name[match(reduced_generic_PPI$TCGA_gene_name, significant_human_protein$TCGA_gene_name)]))
# [1] TRUE
colnames(reduced_generic_PPI)[1:2] = c("Target", "Source")
reduced_generic_PPI$MIST = -1
reduced_generic_PPI$point_col = significant_human_protein$point_col[match(reduced_generic_PPI$Target, significant_human_protein$TCGA_gene_name)]
reduced_generic_PPI$Deviance_hepB = significant_human_protein$Deviance_hepB[match(reduced_generic_PPI$Target, significant_human_protein$TCGA_gene_name)]
reduced_generic_PPI$Deviance_hepB[reduced_generic_PPI$point_col=='decreased_HBV_mut'] = - reduced_generic_PPI$Deviance_hepB[reduced_generic_PPI$point_col=='decreased_HBV_mut']
reduced_generic_PPI$origin = "human_human"
reduced_generic_PPI = reduced_generic_PPI[reduced_generic_PPI$point_col!='non_significant',]
dim(reduced_generic_PPI)
# 85 6
# add edge type - human-viral
viral_and_PPI_nw = cbind(viral_human_interactions, origin="human_viral")
#rename columns
colnames(viral_and_PPI_nw)[1:2] = c("Target", "Source")

aux_red_generic = reduced_generic_PPI[,c(2,1,3,4,5, 6)]
colnames(aux_red_generic) = colnames(reduced_generic_PPI)
all.equal.character(aux_red_generic$Target, as.character(significant_human_protein[match(aux_red_generic$Target, significant_human_protein$TCGA_gene_name),1]))
# [1] TRUE

aux_red_generic = aux_red_generic[aux_red_generic$point_col!='non_significant',]
dim(aux_red_generic)
# [1] 85 6

#create network of viral and human interactions + Reactome ppi interactions that contains significant genes
all.equal.character( colnames(viral_and_PPI_nw), colnames(reduced_generic_PPI))
all.equal.character( colnames(viral_and_PPI_nw), colnames(aux_red_generic))
viral_and_PPI_nw = rbind(viral_and_PPI_nw, reduced_generic_PPI, aux_red_generic)
dim(viral_and_PPI_nw)
# [1] 341    6

# check for duplicated edges
rownames(viral_and_PPI_nw)=NULL

viral_and_PPI_nw = unique(viral_and_PPI_nw)
# [1] 333  6
write.table(viral_and_PPI_nw, "HBV_viral_and_PPI_nw_fdr_10.txt", sep="\t", col.names = T, row.names = F, quote = F)

# no UBC
viral_and_PPI_nw_UBC = viral_and_PPI_nw[-which(viral_and_PPI_nw$Target%in%"UBC" | viral_and_PPI_nw$Source%in%"UBC"), ]
dim(viral_and_PPI_nw_UBC)
# [1] 289   6




HIC_viral_human_PPI_nw = cbind(HIC_viral_human_interactions, origin="human_viral")
colnames(HIC_viral_human_PPI_nw)[1:2] = c("Target", "Source")
dim(HIC_viral_human_PPI_nw)
# [1] 14  6

HIC_viral_and_PPI_nw = rbind(HIC_viral_human_PPI_nw, reduced_generic_PPI, aux_red_generic)
rownames(HIC_viral_and_PPI_nw)=NULL
dim(HIC_viral_and_PPI_nw)
# [1] 172  4
HIC_viral_and_PPI_nw = unique(HIC_viral_and_PPI_nw)

# remove nodes that aren't connected to the viral proteins or other parts of the network
HIC_viral_and_PPI_nw = HIC_viral_and_PPI_nw[!HIC_viral_and_PPI_nw$Target%in%c('FLNA', 'TTN'),]
HIC_viral_and_PPI_nw = HIC_viral_and_PPI_nw[!HIC_viral_and_PPI_nw$Source%in%c('P4HB'), ]
HIC_viral_and_PPI_nw = HIC_viral_and_PPI_nw[!HIC_viral_and_PPI_nw$Source%in%c('PLOD2'), ]
dim(HIC_viral_and_PPI_nw)
# [1] 168  6
write.table(HIC_viral_and_PPI_nw, "HBV_HIC_viral_and_PPI_nw.txt", sep="\t", col.names = T, row.names = F, quote = F)



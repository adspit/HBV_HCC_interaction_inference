# differential analysis on the extended cohort

#load required packages
library(arm)
library(data.table)
library(ggplot2)
library(qvalue)
library(gplots)
library(ggthemes)
# read in the mutation file 
tcga_lich_extended_mutations = read.table("/Users/adriana/Documents/TCGA_LIHC/TCGA_LIHC/data/new_mutations_extended_cohort_lihc_file.txt", sep='\t', stringsAsFactors = F, header = TRUE, row.names = 1) 


dim(tcga_lich_extended_mutations)
# [1]  366 8767


#remove patients that are pos for both HBV and HCV
tcga_lich_extended_mutations = tcga_lich_extended_mutations[-which(tcga_lich_extended_mutations$hepB==1 & tcga_lich_extended_mutations$hepC==1),]


dim(tcga_lich_extended_mutations)
# [1]  359 8767


# read in predefined oncogenes
oncogenes = read.table("/Users/adriana/Documents/TCGA_LIHC/DifferentiallyMutatedGenes/oncogene_tsg.txt", sep = "\t", row.names = 1)
dim(oncogenes)
# [1] 138   1


results_df = data.frame(gene = character(0), mut_type = character(0), Deviance_hepB = numeric(0), pval_hepB = numeric(0), Deviance_hepC = numeric(0), pval_hepC = numeric(0), mutation_rate_overall =numeric(0), rate_p_hepB=numeric(0), rate_p_hepC = numeric(0), rate_p_hepB_neg =numeric(0), rate_p_hepC_neg = numeric(0), rate_n = numeric(0), estimate_coef_hepC=numeric(0), estimate_coef_hepB=numeric(0), stringsAsFactors = FALSE)

for(gene in colnames(tcga_lich_extended_mutations)[1:(ncol(tcga_lich_extended_mutations)-2)])
  
{
  # null model
  fit=bayesglm(tcga_lich_extended_mutations[,gene]~tcga_lich_extended_mutations[,'hepB']+tcga_lich_extended_mutations[,'hepC'],family='binomial')
  # alternative model with hepB only
  fit_hepB =bayesglm(tcga_lich_extended_mutations[,gene]~tcga_lich_extended_mutations[,'hepB'],family='binomial')
  estimate_coef_hepB = fit_hepB$coefficients[[2]]
  # alternative model with hepC only
  fit_hepC = bayesglm(tcga_lich_extended_mutations[,gene]~tcga_lich_extended_mutations[,'hepC'],family='binomial')
  estimate_coef_hepC = fit_hepC$coefficients[[2]]
  
  anovaChisq_hepC_deviance = anova(fit,fit_hepB, test='Chisq')
  Deviance_hepC = anovaChisq_hepC_deviance$Deviance[2]
  pval_hepC = anovaChisq_hepC_deviance[,'Pr(>Chi)'][2]
  
  anovaChisq_hepB_deviance = anova(fit,fit_hepC, test='Chisq')
  Deviance_hepB = anovaChisq_hepB_deviance$Deviance[2]
  pval_hepB = anovaChisq_hepB_deviance[,'Pr(>Chi)'][2]
  
  
  hepB_p_vector=tcga_lich_extended_mutations[which(tcga_lich_extended_mutations$hepB==1),gene]
  hepC_p_vector=tcga_lich_extended_mutations[which(tcga_lich_extended_mutations$hepC==1),gene]
  virus_n_vector=tcga_lich_extended_mutations[tcga_lich_extended_mutations$hepB==0 & tcga_lich_extended_mutations$hepC==0,gene]
  
  hepB_neg = tcga_lich_extended_mutations[which(tcga_lich_extended_mutations$hepB==0),gene] 
  hepC_neg = tcga_lich_extended_mutations[which(tcga_lich_extended_mutations$hepC==0),gene] 
  #
  mutation_rate_overall = length(which(tcga_lich_extended_mutations[,gene]==1))/sum(1-is.na(tcga_lich_extended_mutations[,gene]))
  
  # sum(1-is.na(hepB_p_vector)) represents the number of patients that are positive for hepB and are not NA; analog for hepC
  rate_p_hepB=sum(hepB_p_vector,na.rm=T)/sum(1-is.na(hepB_p_vector))
  rate_p_hepC=sum(hepC_p_vector,na.rm=T)/sum(1-is.na(hepC_p_vector))
  
  rate_p_hepB_neg=sum(hepB_neg,na.rm=T)/sum(1-is.na(hepB_neg))
  rate_p_hepC_neg=sum(hepC_neg,na.rm=T)/sum(1-is.na(hepC_neg))
  
  # rate of mutations for patients that are hepB and hepC negative
  rate_n=sum(virus_n_vector,na.rm=T)/sum(1-is.na(virus_n_vector))
  
  if (is.na(rate_p_hepB) | is.na(rate_p_hepC) | is.na(rate_n)) {next}
  mut_type='unknown'
  if (oncogenes[gene,1] %in% c('Oncogene','Amplification_Oncogene')) {mut_type='Oncogene'}
  if (oncogenes[gene,1] %in% c('TSG','Homozygous_deletion_TSG')) {mut_type='TSG'}
  
  results_df[nrow(results_df)+1, ] <- list(gene, mut_type, Deviance_hepB, pval_hepB, Deviance_hepC, pval_hepC, mutation_rate_overall, rate_p_hepB, rate_p_hepC, rate_p_hepB_neg, rate_p_hepC_neg, rate_n, estimate_coef_hepC, estimate_coef_hepB)
  message(gene,"\t")
}

dim(results_df)
# [1] 8765   12  


results_df$Deviance_hepB = -results_df$Deviance_hepB
results_df$Deviance_hepC = -results_df$Deviance_hepC

significant_nominal_mut_genes_hbv = results_df[results_df$pval_hepB<0.05, ]
dim(significant_nominal_mut_genes)
# [1] 386  14


# keep only those genes for which we have a certain percentage of the cohort mutated
results_thrsh = results_df[which(results_df$gene%in%colnames(tcga_lich_extended_mutations[,colSums(tcga_lich_extended_mutations)>10])),]
dim(results_thrsh)
# [1] 602    8
results_thrsh$qvalue = qvalue(results_thrsh$pval_hepB)$qvalues
results_thrsh = results_thrsh[order(results_thrsh$pval_hepB, results_thrsh$qvalue),]
results_thrsh$Deviance_hepB[results_thrsh$estimate_coef_hepB<0] = - results_thrsh$Deviance_hepB[results_thrsh$estimate_coef_hepB<0]


# keep only necessary column for heatmap
results_thrs_mat = results_thrsh[which(results_thrsh$qvalue<0.2), c(1,3,5,7,8,9,10,11,12,13,14,15)]
# assign gene names as rownames
rownames(results_thrs_mat) = results_thrs_mat$gene

results_df_thrsh_viral_non_viral_mat = results_df_thrsh_viral_non_viral[which(results_df_thrsh_viral_non_viral$qvalue<0.2),]
rownames(results_df_thrsh_viral_non_viral_mat) = results_df_thrsh_viral_non_viral_mat$gene

top_genes = rbind(results_thrs_mat, results_thrsh[results_thrsh$gene%in%setdiff(results_df_thrsh_viral_non_viral_mat$gene, results_thrs_mat$gene),c(1,3,5,7,8,9,10,11,12,13,14,15)])
dim(top_genes)
# [1] 44 12
rownames(top_genes)[41:44] = c( "CPS1",     "ARHGAP44", "GPS2",     "TMEM102")

top_genes$group = "HBV_significant"
top_genes$group[41:44] = "viral_significant"
top_genes$group[which(top_genes$gene%in% intersect(results_df_thrsh_viral_non_viral_mat$gene, results_thrsh$gene[1:40]))] = "HBV_and_viral_significant"
top_genes$color[top_genes$group=="HBV_and_viral_significant"] = '#dc828a'
top_genes$color[top_genes$group=="HBV_significant"] = '#cddc82'
top_genes$color[top_genes$group=="viral_significant"] = '#9182dc'









# log values of mutation frequency 
pdf("heatmap2_all_significant_genes_clustered_log.pdf", height = 6, width = 13, useDingbats = FALSE)

heatmap.2(t(-log(top_genes[,c(9,6,5,4)]+0.01)), trace = "none", dendrogram = "none", density.info = "none",margins = c(9,6), cexRow = 1.2, cexCol = 1.2,  lhei = c(0.5,2), lwid = c(1,3), breaks = seq(0,5,0.01),  col=colorRampPalette(c( "gray96", "#1b54ff"))(500),  key.title = NA, Rowv = FALSE,  cellnote = round(-t(log(top_genes[,c(9,6,5,4)]+0.01)),2), notecol="black",  ColSideColors = top_genes$color, notecex=1.2)
dev.off()



pdf("heatmap2_all_significant_genes_clustered_log_final_version.pdf", height = 6, width = 13, useDingbats = FALSE)

heatmap.2(t(-log(top_genes[,c(9,6,5,4)]+0.01)), trace = "none", dendrogram = "none", density.info = "none",margins = c(9,6), cexRow = 1.2, cexCol = 1.2,  lhei = c(0.5,2), lwid = c(1,3), breaks = seq(0,5,0.01),  col=colorRampPalette(c("#1b54ff", "#b4d9ff", "gray96"))(500),  key.title = NA, Rowv = FALSE,  cellnote = round(t(-log(top_genes[,c(9,6,5,4)]+0.01)),2), notecol="black",  ColSideColors = top_genes$color, notecex=1.2)
dev.off()

aux2_order =  heatmap.2(t(-log(top_genes[,c(9,6,5,4)]+0.01)), trace = "none", dendrogram = "none", density.info = "none",margins = c(9,6), cexRow = 1.2, cexCol = 1.2,  lhei = c(0.5,2), lwid = c(1,3), breaks = seq(0,5,0.01),  col=colorRampPalette(c("#1b54ff", "gray96"))(500),  key.title = NA, Rowv = FALSE,  cellnote = round(t(-log(top_genes[,c(9,6,5,4)]+0.01)),2), notecol="black",  ColSideColors = top_genes$color, notecex=1.2)

rearranged_top_genes = top_genes[aux2_order$colInd,]



pdf("heatmap2_all_significant_genes.pdf", height = 6, width = 13, useDingbats = FALSE)
heatmap.2(t(rearranged_top_genes[,c(9,6,5,4)]), trace = "none", dendrogram = "none", density.info = "none",margins = c(9,6), cexRow = 1.2, cexCol = 1.2,  lhei = c(0.5,2), lwid = c(1,3), breaks = seq(0,0.3,0.001),  col=colorRampPalette(c("gray96","#1b54ff"))(300),  key.title = NA, Rowv = FALSE,Colv = FALSE,  cellnote = round(t(rearranged_top_genes[,c(9,6,5,4)]),2), notecol="black",  ColSideColors = top_genes$color, notecex=1.2)

par(lend = 1)           # square line ends for the color legend
legend(.6,1.20, xpd=TRUE,     # location of the legend on the heatmap plot
       legend = c("HBV and viral significant", "HBV significant", "viral significant"), # category labels
       col = c('#dc828a', "#cddc82", "#9182dc"),  # color key
       lty= 1,             # line style
       lwd = 10,  # line width
)
dev.off()





# matched deviance for the top genes hepB specific genes 
results_df_thrsh_viral_non_viral = results_df_thrsh_viral_non_viral[match(top_genes$gene, results_df_thrsh_viral_non_viral$gene),]
all.equal.character(top_genes$gene, results_df_thrsh_viral_non_viral$gene)
# [1] TRUE
dim(results_df_thrsh_viral_non_viral)
# [1] 44  9
top_genes$Deviance_both = results_df_thrsh_viral_non_viral$Deviance_virus
aux2_order$colInd
# [1] 27 41 20 40 25 26 21 24 22 23 30 28 29 19 36 34 35  8  1  2  3  7 44 42 43 32 18 33 31 39 37 38  6 17 15 16  4  5 14 13 12 11  9 10

aux2_order$colInd = c(aux2_order$colInd[2:length(aux2_order$colInd)], 27)
rearranged_top_genes = top_genes[aux2_order$colInd,]
rearranged_top_genes[,c(15,3,2)]
pdf("deviances_top_genes_significant_genes_mut_direction_3cases_clustered_new.pdf", height = 13, width = 5, useDingbats = FALSE)
heatmap.2(as.matrix(rearranged_top_genes[,c(2,3,15)]), trace = "none", dendrogram = "none", Colv = FALSE, Rowv = FALSE, density.info = "none",lhei = c(0.5,2), lwid = c(1,3), breaks=seq(-14,14,0.1), col=colorRampPalette(c('#36aac1',"gray96","#ff012c"))(280),  key.title = NA, cellnote = round(rearranged_top_genes[,c(2,3,15)],2),  notecol="black", notecex = 1.2)
dev.off()



# mutation rates

r_top_genes = top_genes[rownames(rearranged_top_genes), ]

pdf("heatmap2_all_significant_genes_clustered_040119_2.pdf", height = 13, width = 6, useDingbats = FALSE)

heatmap.2(as.matrix(r_top_genes[,c(4,5,6,9)]*100), trace = "none", dendrogram = "none", density.info = "none",margins = c(9,6), cexRow = 1.2, cexCol = 1.2,  lhei = c(0.5,2), lwid = c(1,3), breaks = seq(0,50,0.1),  col=colorRampPalette(c( "white", "#4f5cff", "#3545ff" ,"#000eb5"))(500),  key.title = NA, Colv = FALSE, Rowv = FALSE,  cellnote = round(r_top_genes[,c(4,5,6,9)]*100,2), notecol="black", notecex=1.2)
dev.off()




cul4 = data.frame(Sample = rownames(tcga_lich_extended_mutations[,c(grep("CUL4", colnames(tcga_lich_extended_mutations)),8766)]), tcga_lich_extended_mutations[,c(grep("CUL4", colnames(tcga_lich_extended_mutations)),8766)]) 
cul4$CUL4A = as.factor(cul4$CUL4A)
cul4$CUL4B = as.factor(cul4$CUL4B)
cls <- hclust(daisy(cul4[,2:3], metric = "gower"))
cul4 = cul4[cls$order,]
cul4[,2] = as.numeric(levels(cul4[,2]))[cul4[,2]]
cul4[,3] = as.numeric(levels(cul4[,3]))[cul4[,3]]



cul4_m = melt(cul4, id.vars = c("Sample", "hepB"))
cul4_m$Sample = factor(cul4_m$Sample, levels = cul4_m$Sample[1:359])
cul4_m$value = factor(cul4_m$value)

library("ggplot2")
library("plyr")
library("reshape2")
library("scales")
library(cluster)

pdf("heatmap_clustered_mutations_cul4.pdf", useDingbats = FALSE, width=13)
ggplot(cul4_m, aes(Sample, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +scale_fill_manual(values=c("lightgrey", "steelblue"))+
  theme_grey(base_size = 7) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_blank())+ coord_fixed(ratio=9)
dev.off()

require(data.table) 
sel = ppp2_m[ppp2_m$value!=0,c("hepB","variable")]
sel = t(xtabs(formula=~hepB + variable, data=sel))
sel = data.frame(HBV_neg=unlist(sel[,1]), HBV_pos=sel[,2])
ppp2 = data.frame(Sample = rownames(tcga_lich_extended_mutations[,c(grep("PPP2", colnames(tcga_lich_extended_mutations)),8766)]), tcga_lich_extended_mutations[,c(grep("PPP2", colnames(tcga_lich_extended_mutations)),8766)]) 
head(ppp2)

ppp2[,2:13] <- lapply(ppp2[,2:13] , factor)


cls <- hclust(daisy(ppp2[,2:13], metric = "gower"))
cls_2 = hclust(daisy(t(ppp2[,2:13]), metric = "gower"))
ppp2 = ppp2[cls$order,]

ppp2_m = melt(ppp2, id.vars = c("Sample", "hepB"))
ppp2_m$Sample = factor(ppp2_m$Sample, levels = ppp2_m$Sample[1:359])
ppp2_m$value = factor(ppp2_m$value)
ppp2_m$hepB = factor(ppp2_m$hepB)
pdf("heatmap_clustered_mutations_ppp2.pdf", useDingbats = FALSE, width=13)
ggplot(ppp2_m, aes(Sample, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +scale_fill_manual(values=c("lightgrey", "steelblue"))+
  theme_tufte(base_family="Helvetica") + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_blank())+ coord_fixed(ratio=9)
dev.off()

pdf("heatmap_clustered_mutations_ppp2_annotation.pdf", useDingbats = FALSE, width=13)
ggplot(ppp2_m, aes(Sample, variable)) + 
  geom_tile(aes(fill = hepB), colour = "white") +scale_fill_manual(values=c("grey", "black"))+
  theme_tufte(base_family="Helvetica") + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_blank())+ coord_fixed(ratio=9)
dev.off()

write.table(sel, "frequency_viral_status_in_PPP2_related_proteins.txt", col.names = TRUE, row.names = TRUE, quote = FALSE)

# pathway enrichment analysis

library(gplots)
library(RCy3)
library(KEGGREST)
library(reshape2)
library(data.table)
library(plyr)
library(Matrix)
library(mgsa)
library(rWikiPathways)
library("RColorBrewer")
library(stringr)
library(GO.db)
library(GenomicFeatures)
library(rtracklayer)

genes <- fread("gencode.v29.chr_patch_hapl_scaff.annotation.gtf", col.names=c("Chromosome", "Source","Type","Start", "End","Score", "Strand", "Phase", "Attributes"), sep="\t", stringsAsFactors=F)

dim(genes)
# [1] 2989084       9

# select only the transcript info
genes <- genes[genes$Type%in%"transcript", ]
dim(genes)
# [1] 226811      9

# split the Attributes column
genes <- with(genes, {cbind(genes[,-9],colsplit(genes$Attributes, pattern=";", names=c("Gene_id", "Transcript_id", "Gene_type", "Gene_name", "Transcript_type", "Transcript_name", "Level", "Transcript_support_level", "Tag", "Havana_gene", "Havana_transcript")))})
genes$Gene_id <- gsub("gene_id ", "",genes$Gene_id)
genes$Transcript_id <- gsub("transcript_id ", "", genes$Transcript_id)
genes$Gene_type <- gsub("gene_type ", "", genes$Gene_type) 
genes$Gene_status <- gsub("gene_status ", "", genes$Gene_status)
genes$Gene_name <- gsub("gene_name ", "", genes$Gene_name)
genes$Transcript_type <- gsub("transcript_type ", "", genes$Transcript_type)
genes$Transcript_name <- gsub("transcript_name ", "", genes$Transcript_name)
genes$Transcript_status <- gsub("transcript_status ", "", genes$Transcript_status)
genes <- genes[,1:14]


unique(genes$Gene_type)
# [1] " \"transcribed_unprocessed_pseudogene\"" " \"unprocessed_pseudogene\""            
# [3] " \"miRNA\""                              " \"lincRNA\""                           
# [5] " \"protein_coding\""                     " \"processed_pseudogene\""              
# [7] " \"antisense\""                          " \"processed_transcript\""              
# [9] " \"snRNA\""                              " \"transcribed_processed_pseudogene\""  
# [11] " \"sense_intronic\""                     " \"misc_RNA\""                          
# [13] " \"TEC\""                                " \"transcribed_unitary_pseudogene\""    
# [15] " \"snoRNA\""                             " \"scaRNA\""                            
# [17] " \"rRNA_pseudogene\""                    " \"unitary_pseudogene\""                
# [19] " \"3prime_overlapping_ncRNA\""           " \"polymorphic_pseudogene\""            
# [21] " \"bidirectional_promoter_lncRNA\""      " \"sense_overlapping\""                 
# [23] " \"pseudogene\""                         " \"rRNA\""                              
# [25] " \"IG_V_pseudogene\""                    " \"scRNA\""                             
# [27] " \"IG_V_gene\""                          " \"IG_C_gene\""                         
# [29] " \"IG_J_gene\""                          " \"sRNA\""                              
# [31] " \"ribozyme\""                           " \"translated_processed_pseudogene\""   
# [33] " \"vaultRNA\""                           " \"TR_C_gene\""                         
# [35] " \"TR_J_gene\""                          " \"TR_V_gene\""                         
# [37] " \"TR_V_pseudogene\""                    " \"TR_D_gene\""                         
# [39] " \"IG_C_pseudogene\""                    " \"non_coding\""                        
# [41] " \"macro_lncRNA\""                       " \"TR_J_pseudogene\""                   
# [43] " \"IG_J_pseudogene\""                    " \"IG_D_gene\""                         
# [45] " \"IG_pseudogene\""                      " \"Mt_tRNA\""                           
# [47] " \"Mt_rRNA\""                           

genes$Gene_type = gsub("\"", "", genes$Gene_type, fixed = TRUE )
genes$Gene_type = gsub(" ", "", genes$Gene_type, fixed = TRUE )

genes <- genes[genes$Gene_type%in%"protein_coding", ]

dim(genes)
# 165792     14

# dealing with isoforms
# when having the same gene id, I select the min start coord and the max end coord
genes <- ddply(genes, .(Chromosome, Source, Type, Strand, Gene_id, Gene_type, Gene_name), summarise,Start=min(Start), End=max(End) )

dim(genes)
# 28957     9

head(genes)
# Chromosome  Source       Type Strand              Gene_id      Gene_type Gene_name  Start    End
# 1       chr1  HAVANA transcript      +  "ENSG00000186092.6" protein_coding   "OR4F5"  65419  71585
# 2       chr1 ENSEMBL transcript      +  "ENSG00000186092.6" protein_coding   "OR4F5"  69055  70108
# 3       chr1  HAVANA transcript      -  "ENSG00000284733.1" protein_coding  "OR4F29" 450703 451697
# 4       chr1  HAVANA transcript      -  "ENSG00000284662.1" protein_coding  "OR4F16" 685679 686673
# 5       chr1  HAVANA transcript      + "ENSG00000187634.11" protein_coding  "SAMD11" 923928 944575
# 6       chr1 ENSEMBL transcript      + "ENSG00000187634.11" protein_coding  "SAMD11" 925741 944581

genes$Gene_id = gsub("\"", "", genes$Gene_id, fixed = TRUE)
genes$Gene_name = gsub("\"", "", genes$Gene_name, fixed = TRUE)
genes$Gene_name = gsub(" ", "", genes$Gene_name, fixed = TRUE)
genes$Gene_id =  gsub("\\..*", "", genes$Gene_id)

sel_ensembl_ids = genes[which(genes$Gene_name%in%top_genes$gene),]
dim(sel_ensembl_ids)

sel_ensembl_ids = ddply(sel_ensembl_ids, .(Chromosome, Type, Strand, Gene_id, Gene_type, Gene_name), summarise,Start=min(Start), End=max(End) )
dim(sel_ensembl_ids)
# [1] 46  8
sel_ensembl_ids = sel_ ensembl_ids[-c(45,46),]
sel_ensembl_ids = sel_ensembl_ids[match(top_genes$gene, sel_ensembl_ids$Gene_name),]
all.equal.character(sel_ensembl_ids$Gene_name, top_genes$gene)
#TRUE


library(mgsa)
wikip=read.delim("wikipathways_data_Homo_sapiens.tab.txt",stringsAsFactors = F)


# match wiki pathways to Ensembl IDs corresponding to our selected genes 
wikip.set=apply(wikip,1,function(x){
  s=unlist(strsplit(x["Ensembl"],","))
  na.omit(match(s,sel_ensembl_ids$Gene_id))
})

ids=sapply(wikip$Url.to.WikiPathways,function(x){y=strsplit(x,":")[[1]];y[length(y)]})

wikip_res_hbv_significant_only = mgsa( which(top_genes$group=="HBV_significant"),wikip.set)
wikip_res_hbv_and_viral_significant = mgsa( which(top_genes$group=="HBV_and_viral_significant"),wikip.set)
wikip_res_viral = mgsa( which(top_genes$group=="viral_significant"),wikip.set)
wikip_res = mgsa( seq(1:nrow(top_genes)),wikip.set)


# prob.m=prob.m[apply(prob.m,1,max,na.rm=T)>.5,apply(prob.m,2,max,na.rm=T)>.9]
 pw.names=wikip.sel$Pathway.Name[match(colnames(prob.m),ids.sel)]

 
 
 data <- data.frame(Gene_set=c("Gene_set1", "Gene_set2", "Gene_set3", "Gene_set4", "Gene_set5"),
                    NES=runif(5, -3, 3),
                    FDR_q.val=runif(5,0,1),
                    No_of_significant_genes=runif(5, 1, 100))
 
 # Plotting
 library(ggplot2) 
 p <- ggplot(data, aes(NES, Gene_set))
 p + geom_point(aes(colour=FDR_q.val, size=No_of_significant_genes)) +
   scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
   geom_vline(xintercept=0, size=0.5, colour="gray50") +
   theme(panel.background=element_rect(fill="gray95", colour="gray95"),
         panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
         panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
         axis.title.y=element_blank()) +
   expand_limits(x=c(-3,3)) +
   scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
   scale_y_discrete(limits=rev(data$Gene_set))
 
 
 
 
 
 
pdf("wikip_pathway_enrichment.pdf",width=20,height=11, useDingbats = FALSE)
heatmap.2(as.matrix(prob.m),trace="none",cexCol = 1.2,cexRow = 1.2,labCol = pw.names,col=colorRampPalette(c("white","#9accff","#005ab3"))(100),margins=c(15,5),density.info="none", srtCol=45, dendrogram = "none", lhei = c(0.5,2.3), lwid = c(1,5), key.title = "Pathway probability", key.par = list(cex=0.8))
dev.off()



prob.m.rearranged = prob.m[match(top_genes$gene, rownames(prob.m)),]
prob.m.rearranged = prob.m.rearranged[aux2_order$colInd,]


pdf("wikip_pathway_enrichment_alt_reordered.pdf",width=20,height=11, useDingbats = FALSE)
heatmap.2(as.matrix(prob.m.rearranged),trace="none",cexCol = 1.2,cexRow = 1.2,labCol = pw.names,col=colorRampPalette(c("white","#9accff","#005ab3"))(100),margins=c(15,5),density.info="none", srtCol=45, dendrogram = "none", lhei = c(0.5,2.3), lwid = c(1,5), key.title = "Pathway probability", key.par = list(cex=0.8), Rowv = FALSE)
dev.off()



# newer version of WikiPathways
downloadPathwayArchive(organism = "Homo sapiens", format ="gmt")



# GO enrichment analysis
entrez=mget(sel_ensembl_ids$Gene_id,org.Hs.egENSEMBL2EG,ifnotfound=NA)
entrez=sapply(entrez,function(x)x[1])
symbol=rep(NA,length(entrez))
symbol[!is.na(entrez)]=unlist(mget(entrez[!is.na(entrez)],org.Hs.egSYMBOL))

go.ann=readGAF("gene_association.goa_human.txt")
go.sets=go.ann@sets[sapply(go.ann@sets,length)>5]

mgsa.res=mclapply(1:nrow(sel_ensembl_ids),function(i){
  symbol.int=na.omit(match(top_genes$gene,go.ann@itemAnnotations$symbol))
  mgsa(symbol.int,go.sets)
},mc.cores=40)


mgsa.res.ann=lapply(mgsa.res,function(x){if(is.null(x))return(NULL);x@setsResults$id=ids;x@setsResults$name=wikip$Pathway.Name;x})
names(mgsa.res.ann)=sel_ensembl_ids$Gene_name
mgsa.res.ann=mgsa.res.ann[!sapply(mgsa.res,function(x)is.null(x))]



sel=results_thrsh$gene[which(results_thrsh$qvalue<0.2)]
sel.int=na.omit(match(sel,go.ann@itemAnnotations$symbol))
go_res = mgsa(sel.int,go.sets)






pdf("go_pathway_enrichment.pdf",width=20,height=11, useDingbats = FALSE)
heatmap.2(as.matrix(prob.m),trace="none",cexCol = 1.2,cexRow = 1.2,labCol = pw.names,col=colorRampPalette(c("white","#9accff","#005ab3"))(100),margins=c(15,5),density.info="none", srtCol=45, dendrogram = "none", lhei = c(0.5,2.3), lwid = c(1,5), key.title = "Pathway probability", key.par = list(cex=0.8))
dev.off()

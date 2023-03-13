#load required packages
library(arm)


# read in TCGA LIHC  data frame 

LIHC_mutations = read.table("/Users/adriana/Documents/TCGA_LIHC/DifferentiallyMutatedGenes/df4glm_LIHC.txt",h=TRUE,sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dim(LIHC_mutations)
# [1]  193 9294

# read in predefined oncogenes
oncogenes = read.table("/Users/adriana/Documents/TCGA_LIHC/DifferentiallyMutatedGenes/oncogene_tsg.txt", sep = "\t", row.names = 1)
dim(oncogenes)
# [1] 138   1

results_df = data.frame(gene = character(0), mut_type = character(0), Deviance_hepB = numeric(0), pval_hepB = numeric(0), Deviance_hepC = numeric(0), pval_hepC = numeric(0), rate_p_hepB=numeric(0), rate_p_hepC = numeric(0), rate_n = numeric(0), mut_rate = numeric(0), estimate_coef_hepC=numeric(0), estimate_coef_hepB=numeric(0), stringsAsFactors = FALSE)

for(gene in colnames(LIHC_mutations)[1:(ncol(LIHC_mutations)-2)])
  
{
  fit=bayesglm(LIHC_mutations[,gene]~LIHC_mutations[,'hepB']+LIHC_mutations[,'hepC'],family='binomial')
  
  
  # 
  # summary(fit)
  # 
  # Call:
  #   glm(formula = LIHC_mutations[, gene] ~ LIHC_mutations[, "hepB"] + 
  #         LIHC_mutations[, "hepC"], family = "binomial")
  # 
  # Deviance Residuals: 
  #   Min       1Q   Median       3Q      Max  
  # -0.3300  -0.2345  -0.2137  -0.2137   2.7537  
  # 
  # Coefficients:
  #   Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)               -3.7686     0.5592  -6.739  1.6e-11 ***
  #   LIHC_mutations[, "hepB"]   0.6963     1.1434   0.609    0.543    
  # LIHC_mutations[, "hepC"]   0.1890     1.1369   0.166    0.868    
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # (Dispersion parameter for binomial family taken to be 1)
  # 
  # Null deviance: 46.402  on 192  degrees of freedom
  # Residual deviance: 46.058  on 190  degrees of freedom
  # AIC: 52.058
  # 
  # Number of Fisher Scoring iterations: 6
  
  fit_hepB =bayesglm(LIHC_mutations[,gene]~LIHC_mutations[,'hepB'],family='binomial')
  estimate_coef_hepB = fit_hepB$coefficients[[2]]
  fit_hepC = bayesglm(LIHC_mutations[,gene]~LIHC_mutations[,'hepC'],family='binomial')
  estimate_coef_hepC = fit_hepC$coefficients[[2]]
  
  anovaChisq_hepC_deviance = anova(fit,fit_hepB, test='Chisq')
  # print(anovaChisq_hepC_deviance)
  # Analysis of Deviance Table
  # 
  # Model 1: LIHC_mutations[, gene] ~ LIHC_mutations[, "hepB"] + LIHC_mutations[,"hepC"]
  # Model 2: LIHC_mutations[, gene] ~ LIHC_mutations[, "hepB"]
  # Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
  # 1       190     46.058                      
  # 2       191     46.085 -1 -0.026657   0.8703
  # deviance between the model in which we use both hepC and hepB status and the model where we use only the hepB as a covariate
  Deviance_hepC = anovaChisq_hepC_deviance$Deviance[2]
  pval_hepC = anovaChisq_hepC_deviance[,'Pr(>Chi)'][2]
  
  
  anovaChisq_hepB_deviance = anova(fit,fit_hepC, test='Chisq')
  # print(anovaChisq_hepB_deviance)
  # Analysis of Deviance Table
  # 
  # Model 1: LIHC_mutations[, gene] ~ LIHC_mutations[, "hepB"] + LIHC_mutations[, "hepC"]
  # Model 2: LIHC_mutations[, gene] ~ LIHC_mutations[, "hepC"]
  # Resid. Df Resid. Dev Df Deviance Pr(>Chi)
  # 1       190     46.058                     
  # 2       191     46.382 -1  -0.3245   0.5689
  # deviance between the model in which we use both hepC and hepB status and the model where we use only the hepC as a covariate
  Deviance_hepB = anovaChisq_hepB_deviance$Deviance[2]
  pval_hepB = anovaChisq_hepB_deviance[,'Pr(>Chi)'][2]
  
  
  hepB_p_vector=LIHC_mutations[which(LIHC_mutations$hepB==1),gene]
  hepC_p_vector=LIHC_mutations[which(LIHC_mutations$hepC==1),gene]
  virus_n_vector=LIHC_mutations[LIHC_mutations$hepB==0 & LIHC_mutations$hepC==0,gene]
  
  # sum(1-is.na(hepB_p_vector)) represents the number of patients that are positive for hepB and are not NA; analog for hepC
  rate_p_hepB=sum(hepB_p_vector,na.rm=T)/sum(1-is.na(hepB_p_vector))
  rate_p_hepC=sum(hepC_p_vector,na.rm=T)/sum(1-is.na(hepC_p_vector))
  
  # rate of mutations for patients that are hepB and hepC negative
  rate_n=sum(virus_n_vector,na.rm=T)/sum(1-is.na(virus_n_vector))
  
  if (is.na(rate_p_hepB) | is.na(rate_p_hepC) | is.na(rate_n)) {next}
  
  mut_rate=sum(LIHC_mutations[,gene],na.rm=T)/sum(1-is.na(LIHC_mutations[,gene]))
  mut_type='loss'
  if (oncogenes[gene,1] %in% c('Oncogene','Amplification_Oncogene')) {mut_type='gain'}
  
  results_df[nrow(results_df)+1, ] <- list(gene, mut_type, -Deviance_hepB, pval_hepB, -Deviance_hepC, pval_hepC, rate_p_hepB, rate_p_hepC, rate_n, mut_rate, estimate_coef_hepC, estimate_coef_hepB)
  message(gene,"\t")
}

results_df = results_df[order(results_df$pval_hepC, results_df$Deviance_hepC),]

write.table(results_df, "diff_mutated_genes_HCV_bayes_lm_TCGA_LIHC_extended.txt", col.names = TRUE, sep='\t', row.names = FALSE, quote = FALSE)




top_genes_hepC = results_df[order(results_df$pval_hepC, results_df$pval_hepB),]
top_genes_hepC = top_genes_hepC[top_genes_hepC$pval_hepC < 0.05,]
dim(top_genes_hepC)
# [1] 248  12

top_genes_hepB = results_df[order(results_df$pval_hepB, results_df$pval_hepC),]
top_genes_hepB = top_genes_hepB[top_genes_hepB$pval_hepB < 0.05,]
dim(top_genes_hepB)
# [1] 311 12



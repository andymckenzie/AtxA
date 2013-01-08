#script is distributed under the mit license, http://www.opensource.org/licenses/MIT, with no warranty of any kind
#given a set of genes with replicates for both environment and genotype
#calls differential expression based solely on environment and genotype, via likelihood ratio test
#calls differential expression, based on environment by genotype interaction, via likelihood ratio test
#calls fold change of the above comparisons
#generates a png of volcano plot to summarize the joint distribution of the above calls 
#primary dependency is the R package bioconductor

#load packages
library ("DESeq") 
library (ggplot2) 

#specify count file, where rows = genes, columns = conditions, and cells = raw counts 
countfile = "chrom_counts.txt" 

#specify normalized count file, where rows = genes, columns = conditions, and cells = normalized counts
normcountfile = "norm_chrom_counts.txt"

#specify expression quantile (across all conditions; total row sum count) below which genes are not included
expression_cutoff = 0.3

#specify p-value threshold
p_value_threshold = 0.01

#specify upregulation of main effect necessary to be called as upregulated
fold_change_main_cutoff = 2

#specify upregulation of interaction effect (environment-by-genotype pairing) necessary to be called as upregulated
fold_change_interaction_cutoff = 4

#read the count file into a table
atxaCountTable = read.table (countfile, header=TRUE, row.names=1)

#read the normalized count file into a table
normCountTable = read.table (normcountfile, header=TRUE, row.names=1)

#specify environment and genotype assignments for each sample
atxaDesign = data.frame(row.names = colnames(atxaCountTable), enviro = c("CO2", "CO2", "air", "air", "air", "air", "CO2", "CO2"), genotype = c("Ames35", "Ames35", "dAtxA", "dAtxA", "Ames35", "Ames35", "dAtxA", "dAtxA"))

#specify design matrix of experiment
cdsFull = newCountDataSet (atxaCountTable, atxaDesign)

#estimate the size factors for the design matrix
cdsFull = estimateSizeFactors (cdsFull)

#estimate the dispersion parameters for the design matrix
cdsFull = estimateDispersions (cdsFull) 

#filter out genes with expression below the quantile cutoff from the normalized count table
rowsum_norm = rowSums (normCountTable) 
use_filt = (rowsum_norm > quantile(rowsum_norm, expression_cutoff))
cdsFilt = cdsFull [use_filt,]
normCountTable_Filt = normCountTable [use_filt,]  

#generate fits of the generalized linear model for each gene, genotype only
fit0a = fitNbinomGLMs(cdsFilt, count ~ genotype, glmControl = list(maxit = 100))

#generate fits of the generalized linear model for each gene, environment only
fit0b = fitNbinomGLMs(cdsFilt, count ~ enviro, glmControl = list(maxit = 100))

#generate fits of the generalized linear model for each gene, excluding the interaction
fit1 = fitNbinomGLMs(cdsFilt, count ~ genotype + enviro, glmControl = list(maxit = 100))

#generate fits of the generalized linear model for each gene, including the interaction
fit2 = fitNbinomGLMs(cdsFilt, count ~ genotype + enviro + genotype:enviro, glmControl = list(maxit = 100))

#compare fit 1 to fit0b to test whether the genotype has an effect
pvalsGLM_genotype = nbinomGLMTest(fit1, fit0b)

#compare fit 1 to fit0a to test whether the environment has an effect
pvalsGLM_environment = nbinomGLMTest(fit1, fit0a)

#compare fit2 to fit1 to test the effect of interaction
pvalsGLM_interaction = nbinomGLMTest(fit2, fit1)

#control false discovery rates with benjamini-hochberg
padjGLM_genotype = p.adjust(pvalsGLM_genotype, method="BH") 
padjGLM_environment = p.adjust(pvalsGLM_environment, method="BH") 
padjGLM_interaction = p.adjust(pvalsGLM_interaction, method="BH") 

#calculate the mean normalized expression value for each environment and genotype pairing
Ames35_CO2_mean = apply(cbind(normCountTable_Filt$Ames35_CO2_1, normCountTable_Filt$Ames35_CO2_2), 1, mean) 
Ames35_air_mean = apply(cbind(normCountTable_Filt$Ames35_Air_1, normCountTable_Filt$Ames35_Air_2), 1, mean) 
dAtxA_CO2_mean = apply(cbind(normCountTable_Filt$Ames35delAtxA_CO2_1, normCountTable_Filt$Ames35delAtxA_CO2_2), 1, mean) 
dAtxA_air_mean = apply(cbind(normCountTable_Filt$Ames35delAtxA_Air_1, normCountTable_Filt$Ames35delAtxA_Air_2), 1, mean)

#calculate fold changes in CO2 vs air, all cells
fold_change_CO2_air = apply(cbind(Ames35_CO2_mean, dAtxA_CO2_mean), 1, mean) / apply(cbind(Ames35_air_mean, dAtxA_air_mean), 1, mean)

#calculate fold changes in WT vs dAtxA, all cells
fold_change_WT_dAtxA = apply(cbind(Ames35_CO2_mean, Ames35_air_mean), 1, mean) / apply(cbind(dAtxA_CO2_mean, dAtxA_air_mean), 1, mean)

#calculate fold changes in Ames35_CO2 vs dAtxa_air
fold_change_WTCO2_dAtxAair = Ames35_CO2_mean / dAtxA_air_mean

#calculate fold changes in Ames35_air vs dAtxA_CO2
fold_change_WTair_dAtxACO2 = Ames35_air_mean / dAtxA_CO2_mean

#combine the p-value and fold-change calls into one table
table_genotype = cbind(row.names(normCountTable_Filt), padjGLM_genotype, fold_change_WT_dAtxA) 
table_environment = cbind(row.names(normCountTable_Filt), padjGLM_environment, fold_change_CO2_air) 
table_WTCO2_dAtxAair = cbind(row.names(normCountTable_Filt), padjGLM_interaction, fold_change_WTCO2_dAtxAair) 
table_WTair_dAtxACO2 = cbind(row.names(normCountTable_Filt), padjGLM_interaction, fold_change_WTair_dAtxACO2)

#parse the main effect tables into groups of genes that are significantly different and pass the expression filters
signif_upreg_WT = subset(table_genotype, padjGLM_genotype < p_value_threshold & fold_change_WT_dAtxA > fold_change_main_cutoff)
signif_upreg_dAtxA = subset(table_genotype, padjGLM_genotype < p_value_threshold & fold_change_WT_dAtxA < (1/fold_change_main_cutoff))
signif_upreg_CO2 = subset(table_environment, padjGLM_environment < p_value_threshold & fold_change_CO2_air > fold_change_main_cutoff)
signif_upreg_air = subset(table_environment, padjGLM_environment < p_value_threshold & fold_change_CO2_air < (1/fold_change_main_cutoff))

#parse the interaction effect tables into groups of genes that are significantly different and pass the expression filters
signif_upreg_WTCO2 = subset(table_WTCO2_dAtxAair, padjGLM_interaction < p_value_threshold & fold_change_WTCO2_dAtxAair > fold_change_interaction_cutoff)
signif_upreg_dAtxAair = subset(table_WTCO2_dAtxAair, padjGLM_interaction < p_value_threshold & fold_change_WTCO2_dAtxAair < (1/fold_change_interaction_cutoff))
signif_upreg_WTair = subset(table_WTair_dAtxACO2, padjGLM_interaction < p_value_threshold & fold_change_WTair_dAtxACO2 > fold_change_interaction_cutoff)
signif_upreg_dAtxACO2 = subset(table_WTair_dAtxACO2, padjGLM_interaction < p_value_threshold & fold_change_WTair_dAtxACO2 < (1/fold_change_interaction_cutoff))

#write the above comparisons to table in external file
write.table(signif_upreg_WT, file = "signif_upreg_WT.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_dAtxA, file = "signif_upreg_dAtxA.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_CO2, file = "signif_upreg_CO2.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_air, file = "signif_upreg_air.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_WTCO2, file = "signif_upreg_WTCO2.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_dAtxAair, file = "signif_upreg_dAtxAair.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_WTair, file = "signif_upreg_WTair.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(signif_upreg_dAtxACO2, file = "signif_upreg_dAtxACO2.txt", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE) 

#highlight genes that are in the "extended pathogenicity island" for the volcano plots
normCountTable_Filt$threshold = as.factor(row.names(normCountTable_Filt) >= "GBAA_pXO1_0122" & row.names(normCountTable_Filt) <= "GBAA_pXO1_0194" | row.names(normCountTable_Filt) == "sRNA-1" | row.names(normCountTable_Filt) == "sRNA-2")

#construct the modified volcano plot 
g = ggplot(data=normCountTable_Filt, aes(x=log2(fold_change_WTCO2_dAtxAair), y=-log10(padjGLM_interaction), colour=threshold)) 
png(filename = "chrom_counts.png", width = 1000, height = 1000, units = 'px')
print(g +
  geom_point(alpha=0.95, size=4.5) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) + 
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("Log2 Fold Change") + ylab("-Log10 P-Value") +
  geom_hline(aes(yintercept=2), linetype = 2, colour = "#008080", size = 3.5) + 
  theme_bw() + opts(legend.position="none") 
)
dev.off()

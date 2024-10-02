if (!requireNamespace("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
if (!requireNamespace("viridis", quietly = TRUE)){install.packages("viridis")}
if (!requireNamespace("wacolors", quietly = TRUE)){install.packages("wacolors")}

library("wacolors")
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("viridis")
library("plyr")

#set wd and create the necessary ones for results
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/plots/gene_counts/",sep = "")
dir.create(resultsPath)
# now get a summary of gene expression separated by groups, AD vs Non-AD, are 
#there any genes expressed in most samples in one group but not in the other?
#let's see
counts_binary <- counts_matrix
counts_binary[counts_binary>0] <- 1
exp_in_samples_AD   <- data.frame(n_exp = rowSums(counts_binary[,which(annotation$AD=='AD')],na.rm = TRUE))
exp_in_samples_NoAD <- data.frame(n_exp = rowSums(counts_binary[,which(annotation$AD=='No_AD')],na.rm = TRUE))

file_name <- paste(resultsPath,"/ROSMAP_number_of_expressions_per_gene_NoAD.pdf",sep="")
p <- ggplot(exp_in_samples_NoAD, aes(x=n_exp)) + 
     geom_histogram(binwidth = 1) + scale_fill_wa_d(wacolors$volcano) + 
     theme_minimal() + xlab("Number of ocurrences")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

file_name <- paste(resultsPath,"/ROSMAP_number_of_expressions_per_gene_AD.pdf",sep="")
p <- ggplot(exp_in_samples_AD, aes(x=n_exp)) + 
     geom_histogram(binwidth = 1) + scale_fill_wa_d(wacolors$volcano) + 
     theme_minimal() + xlab("Number of ocurrences")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
#The number of ocurrences per gene follows a U-shaped distribution for both AD and 
#non-AD samples, let's identify those genes that are never expressed
neverExpressed <- which(exp_in_samples_AD==0 & exp_in_samples_NoAD==0)
noExpGenes <- as.data.frame(rownames(counts_matrix[neverExpressed,]))
#save the list of never expressed genes (maybe for GSEA for control, basically non
#neuronal gene sets should pop-up here) and also remove them from the dataset for 
#further analysis
counts_matrix <- counts_matrix[-neverExpressed,]
exp_in_samples_AD <- exp_in_samples_AD[-neverExpressed]
exp_in_samples_NoAD <- exp_in_samples_NoAD[-neverExpressed]
write_delim(noExpGenes, file = 'data/genes_never_expressed.txt',delim = '\t', na='NA')
#Keep genes expressed in at least 25% of the samples
low_occurrence_genes_noAD <- which(exp_in_samples_NoAD<0.25*ncol(counts_matrix))
low_occurrence_genes_AD <- which(exp_in_samples_AD<0.25*ncol(counts_matrix))
#check for genes that are lowly ocurring in both subsets
low_occur <- intersect(low_occurrence_genes_noAD,low_occurrence_genes_AD)
#Are there highly occurring genes in AD that are low occuring in non-AD?
high_occur_AD <- which(exp_in_samples_AD>=0.75*ncol(counts_matrix))
opposite_occur_AD <- intersect(high_occur_AD,low_occurrence_genes_noAD)
#there are no such cases, test for the opposite then
high_occur_NoAD <- which(exp_in_samples_NoAD>=0.75*ncol(counts_matrix))
opposite_occur_NoAD <- intersect(high_occur_NoAD,low_occurrence_genes_AD)
#Same, there seems to be no opposite occurence of gene expression among AD vs Non-AD samples
#Now let's be more granular and analyse this sample by sample (AD)
AD_idxs <- which(annotation$AD)
AD_sample_exclusive_genes <- data.frame(id = numeric(length(AD_idxs)), exc_genes = numeric(length(AD_idxs)), n_exc_genes = numeric(length(AD_idxs)))
AD_sample_exclusive_genes$exc_genes   <- 0
AD_sample_exclusive_genes$n_exc_genes <- 0
counter <- 1
for (i in AD_idxs){
  AD_sample_i <- counts_matrix[,i]
  difference  <- which(AD_sample_i>0 & exp_in_samples_NoAD==0)
  n_genes     <- length(difference)
  AD_sample_exclusive_genes$id[counter] <- i
  if (n_genes>0){
    AD_sample_exclusive_genes$exc_genes[counter]   <- paste(difference,collapse=",")
    AD_sample_exclusive_genes$n_exc_genes[counter] <- n_genes
  }
  counter <- counter+1
}
#visualize this on a histogram
#pdf(file=file_name,width=4,height=4)
p <- ggplot(AD_sample_exclusive_genes, aes(x=n_exc_genes)) + 
     geom_histogram(binwidth = 1) + scale_fill_wa_d(wacolors$volcano) + theme_minimal()
file_name <- paste(resultsPath,"/ROSMAP_AD_excl_genes_perSample.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

#the majority of AD samples express less than 5 genes that are not expressed in the non-AD samples
#let's identify these genes, and check their frequency of expression in AD samples
AD_gene_idx_vector <- as.numeric(unlist(strsplit(paste(AD_sample_exclusive_genes$exc_genes,collapse=","), ',')))
AD_gene_idx_vector <- AD_gene_idx_vector[AD_gene_idx_vector>0]
#AD_gene_idx_vector <- data.frame(genes = AD_gene_idx_vector)
AD_gene_excl_freq <- data.frame(id = rownames(counts_matrix)[unique(AD_gene_idx_vector)],
                                index = unique(AD_gene_idx_vector),
                                freq  = rep(0,length(unique(AD_gene_idx_vector))))
for (i in 1:length(AD_gene_excl_freq$index)){
  AD_gene_excl_freq$freq[i] <- length(which(AD_gene_idx_vector == AD_gene_excl_freq$index[i]))
}
p <- ggplot(AD_gene_excl_freq, aes(x=freq)) + 
  geom_histogram(binwidth = 1) + scale_fill_wa_d(wacolors$volcano) + theme_minimal()
file_name <- paste(resultsPath,"/ROSMAP_AD_excl_genes_frequency.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
#get those AD-specific genes that are expressed in more than 1 sample
AD_exclusive_genes_repeated <- AD_gene_excl_freq[AD_gene_excl_freq$freq>1,]
write_delim(AD_exclusive_genes_repeated, file = 'results/AD_exclusiveGenes_repeated.txt',delim = '\t', na='NA')
write_delim(data.frame(id = rownames(counts_matrix)), file = 'results/ROSMAP_expressed_genes.txt',delim = '\t', na='NA')

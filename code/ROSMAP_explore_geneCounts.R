if (!requireNamespace("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
if (!requireNamespace("viridis", quietly = TRUE)){install.packages("viridis")}
if (!requireNamespace("wacolors", quietly = TRUE)){install.packages("wacolors")}

library("wacolors")
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("viridis")
#set wd and create the necessary ones for results
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/plots/gene_counts/",sep = "")
dir.create(resultsPath)
#load data
gene_ids <- read_delim(file = 'data/ROSMAP_annotated_samples_geneIDs.txt',delim = '\t', na='NA')
counts_matrix <- read_delim(file = 'data/ROSMAP_annotated_samples_counts.txt',delim = '\t', na='NA')
annotation    <- read_delim(file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')
counts_matrix <- as.data.frame(counts_matrix)
annotation$braaksc %<>% as.numeric
rownames(counts_matrix) <- t(gene_ids)
# analyze distributions per batch
#let's see the batches at disposal
batches <- as.factor(annotation$Batch)
print(batches)
#Samples belonging to multiple batches were found, let's delete them from the study
idx2rmv       <- which(annotation$Batch == "0, 6, 7")
annotation    <- annotation[-idx2rmv,]
counts_matrix <- counts_matrix[,-idx2rmv]
df_counts     <- as.data.frame(t(counts_matrix))
df_counts$batch <- annotation$Batch
df_counts$AD    <- annotation$AD


counts2plot<- data.frame(counts = c(t(df_counts[,1:(ncol(df_counts)-2)])),batch = rep(df_counts[,ncol(df_counts)-1],nrow(counts_matrix)),AD = rep(df_counts[,ncol(df_counts)],nrow(counts_matrix)))#, AD = c(df_counts[,ncol(df_counts)]))
counts2plot$batch <- as.factor(counts2plot$batch)
counts2plot$counts <- as.numeric(counts2plot$counts)
counts2plot$AD <- as.factor(counts2plot$AD)
counts2plot$cpm <- counts2plot$counts/1E6#log10(counts2plot$counts+1)
counts2plot$lcpm <- log2(counts2plot$cpm)
#counts2plot$counts[counts2plot$counts>2] <- 2
p <- ggplot(counts2plot, aes(x=batch, y=lcpm, fill=AD)) +
  geom_violin(trim=TRUE) + theme_minimal() + scale_fill_wa_d(wacolors$volcano)   #+
  #stat_summary(fun.data="mean_sdl", mult=1,geom="crossbar", width=0.2) +
  #stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="red")
file_name <- paste(resultsPath,"/ROSMAP_counts_perBatch.pdf",sep="")
ggsave(file_name, p, width = 12, height = 10, units = "cm",dpi = 400)
#Distributions of gene counts for batch 0 look a bit different to the rest, let's test 
#using Kolmogorov-Smirnov test
test_matrix <- matrix(1, 9, 9) 
for (i in 0:8){
  dist_i <- counts2plot$cpm[counts2plot$batch==i]
  for (j in 0:8){
    dist_j   <- counts2plot$cpm[counts2plot$batch==j]
    #if(j<i){
      KSresult <- ks.test(dist_i,dist_j, alternative = "two.sided")
      test_matrix[i+1,j+1] <- KSresult$p.value
    #}
  }
}
#significant differences were found among almost all pairwise comparisons, not 
#possible to draw conclusions from this

# now get a summary of gene expression separated by groups, AD vs Non-AD, are 
#there any genes expressed in most samples in one group but not in the other?
#let's see
#AD_counts    <- counts_matrix[,which(annotation$AD==TRUE)]
#nonAD_counts <- counts_matrix[,which(annotation$AD==FALSE)]
counts_binary <- counts_matrix
counts_binary[counts_binary>0] <- 1
exp_in_samples_AD <- rowSums(counts_binary[,which(annotation$AD==TRUE)],na.rm = TRUE)
exp_in_samples_NoAD <- rowSums(counts_binary[,which(annotation$AD==FALSE)],na.rm = TRUE)

file_name <- paste(resultsPath,"/ROSMAP_number_of_expressions_per_gene_NoAD.pdf",sep="")
pdf(file=file_name,width=4,height=4)
hist(exp_in_samples_NoAD, main = "Expressed in non-AD samples",
          xlim = c(0,length(which(annotation$AD==FALSE))),xlab = "Number of ocurrences")
dev.off()

file_name <- paste(resultsPath,"/ROSMAP_number_of_expressions_per_gene_AD.pdf",sep="")
pdf(file=file_name,width=4,height=4)
hist(exp_in_samples_AD, main = "Expressed in AD samples",
     xlim = c(0,length(which(annotation$AD==TRUE))),xlab = "Number of ocurrences")
dev.off()

#The number of ocurrences per gene follows a U-shaped distribution for both AD and 
#non-AD samples, let's identify those genes that are never expressed
neverExpressed <- which(exp_in_samples_AD==0 & exp_in_samples_NoAD==0)
noExpGenes <- as.data.frame(rownames(counts_matrix[neverExpressed,]))
#save the list of never expressed genes (maybe for GSEA for control, basically non
#neuronal gene sets should pop-up here) and also remove them from the dataset for 
#further analysis
temp <- counts_matrix[-noExpGenes]
write_delim(noExpGenes, file = 'data/genes_never_expressed.txt',delim = '\t', na='NA')
#Keep genes expressed in at least 20% of the samples
low_occurrence_genes_noAD <- which(exp_in_samples_NoAD<0.2*ncol(counts_matrix))
low_occurrence_genes_AD <- which(exp_in_samples_AD<0.2*ncol(counts_matrix))
low_occur <- intersect(low_occurrence_genes_noAD,low_occurrence_genes_AD)
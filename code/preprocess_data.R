if (!require("BiocManager", quietly = TRUE))
  {install.packages("BiocManager")}
#BiocManager::install("sva")
if (!require("sva", quietly = TRUE))
{BiocManager::install("sva")}
if (!requireNamespace("ggplot2", quietly = TRUE)){install.packages("ggfortify")}
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
library('sva')
library("ggfortify")
#set wd and create the necessary ones for results
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/DE_analysis/",sep = "")
dir.create(resultsPath)
#load data
gene_ids <- read_delim(file = 'data/ROSMAP_annotated_samples_geneIDs.txt',delim = '\t', na='NA')
genesNoVersion <-  gsub(pattern = '.[0−9]*$',x =gene_ids$gene_ids, replacement = '')
genesNoVersion <- unique(genesNoVersion)
counts_matrix <- read_delim(file = 'data/ROSMAP_annotated_samples_counts.txt',delim = '\t', na='NA')
annotation    <- read_delim(file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')
counts_matrix <- as.data.frame(counts_matrix)
rownames(counts_matrix) <- t(gene_ids)
# analyze distributions per batch
#let's see the batches
batches <- as.factor(annotation$Batch)
#
#Samples belonging to multiple batches were found, let's delete them from the study
idx2rmv       <- which(annotation$Batch == "0, 6, 7")
annotation    <- annotation[-idx2rmv,]
counts_matrix <- counts_matrix[,-idx2rmv]
#delete thos gene entries with 0 counts across all samples
idx2rmv       <- which(rowSums(counts_matrix)==0)
counts_matrix <- counts_matrix[-idx2rmv,]
df_counts     <- as.data.frame(t(counts_matrix))
df_counts$batch <- annotation$Batch
annotation$AD[is.na(annotation$AD)] <- 'other'
df_counts$AD <- annotation$AD
#focus on 
# p <- ggplot(annotation,aes(x = Batch,fill = AD)) + theme_bw() +
#      geom_bar(position = "stack") + scale_fill_wa_d(wacolors$volcano) + 
#      xlab("Sequencing batches") + ylab("Number of samples")
# file_name <- paste(resultsPath,"/ROSMAP_batch_distribution.pdf",sep="")
# ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
# 
# counts2plot<- data.frame(counts = c(t(df_counts[,1:(ncol(df_counts)-2)])),batch = rep(df_counts[,ncol(df_counts)-1],nrow(counts_matrix)),AD = rep(df_counts[,ncol(df_counts)],nrow(counts_matrix)))#, AD = c(df_counts[,ncol(df_counts)]))
# counts2plot$batch <- as.factor(counts2plot$batch)
# counts2plot$counts <- as.numeric(counts2plot$counts)
# counts2plot$AD <- as.factor(counts2plot$AD)
# counts2plot$cpm <- counts2plot$counts/(1E6)#log10(counts2plot$counts+1)
# counts2plot$lcpm <- log2(counts2plot$cpm)
# #counts2plot$counts[counts2plot$counts>2] <- 2
# p <- ggplot(counts2plot, aes(x=batch, y=lcpm, fill=AD)) +
#   geom_violin(trim=TRUE) + theme_minimal() + scale_fill_wa_d(wacolors$volcano) #+
# #stat_summary(fun.data="mean_sdl", mult=1,geom="crossbar", width=0.2) +
# #stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="red")
# file_name <- paste(resultsPath,"/ROSMAP_counts_perBatch.pdf",sep="")
# ggsave(file_name, p, width = 12, height = 10, units = "cm",dpi = 400)
# 
t_c_matrix   <- t(counts_matrix)
t_c_matrix <- as.data.frame(t(counts_matrix))
t_c_matrix$batch <- annotation$Batch
#pca_object   <- prcomp(t_c_matrix, center = TRUE, scale. = TRUE)
#autoplot(pca_object, data= t_c_matrix, colour = 'batch')

#remove batch 7 and 0
toRmv      <- which(annotation$Batch==7 | annotation$Batch==0)
annotation <- annotation[-toRmv,]
AD_col     <- annotation$AD
batch_col  <- annotation$Batch
counts_matrix_f <- counts_matrix[,-toRmv]
idx2rmv         <- which(rowSums(counts_matrix_f)==0)
counts_matrix_f <- counts_matrix_f[-idx2rmv,]
t_c_matrix      <- as.numeric(t(counts_matrix_f))
t_c_matrix      <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- batch_col
#pca_object      <- prcomp(t_c_matrix, center = TRUE, scale. = TRUE)
#autoplot(pca_object, data= t_c_matrix, colour = 'batch')
#cpm_matrix_f <- counts_matrix_f/(colSums(counts_matrix_f)/1E6)
sample_ID <- AD_col
sample_bt <- batch_col
toRemove <- which(AD_col=='other')
annotation <-annotation[-toRemove,]
sample_ID <- sample_ID[-toRemove]
sample_bt <- sample_bt[-toRemove]

counts_matrix.processed <- counts_matrix_f[,-toRemove]
sample_ID <- paste(sample_ID,1:length(sample_ID),sep='_')
sample_ID[which(grepl('No_AD_',sample_ID))] <- 'NCI'
sample_ID <- as.data.frame(sample_ID)
orderPos <- order(sample_ID$sample_ID)
annotation <- annotation[orderPos,]
sample_ID <- sample_ID[orderPos,]
sample_bt <- sample_bt[orderPos]
counts_matrix.processed <- counts_matrix.processed[,orderPos]
#remove non-expressed genes
idx2rmv <- which(rowSums(counts_matrix.processed)==0)
counts_matrix.processed <- counts_matrix.processed[-idx2rmv,]
x <- factor(sample_ID)
#avoid NA and negative values
counts_matrix.processed[counts_matrix.processed<0] <- 0
counts_matrix.processed[is.na(counts_matrix.processed)] <- 0
#PCA
t_c_matrix       <- t(counts_matrix.processed)
#t_c_matrix$AD <- 
pca_object       <- prcomp(t_c_matrix, center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(sample_bt)
autoplot(pca_object, data= t_c_matrix, colour = 'batch')
#remove two oulier samples (PC2) from batch 5. (145_120419)
x    <- pca_object$x
pos  <- which(x[,1]==max(x[,1]))
pos2 <- which(x[,2]==min(x[,2]))
counts_matrix.processed <- counts_matrix.processed[,-pos]
annotation <- annotation[-pos,]
sample_ID  <- sample_ID[-pos]
sample_bt  <- sample_bt[-pos]
idx2rmv    <- which(rowSums(counts_matrix.processed)==0)
counts_matrix.processed <- counts_matrix.processed[-idx2rmv,]
t_c_matrix       <- t(counts_matrix.processed)
pca_object       <- prcomp(as.matrix(t_c_matrix), center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(sample_bt)
autoplot(pca_object, data= t_c_matrix, colour = 'batch')
#remove last outlier. 122_120418
x <- pca_object$x
pos2 <- which(x[,2]==min(x[,2]))
counts_matrix.processed <- counts_matrix.processed[,-pos2]
annotation <- annotation[-pos2,]
sample_ID <- sample_ID[-pos2]
sample_bt <- sample_bt[-pos2]
idx2rmv         <- which(rowSums(counts_matrix.processed)==0)
counts_matrix.processed <- counts_matrix.processed[-idx2rmv,]
t_c_matrix     <- t(counts_matrix.processed)
pca_object      <- prcomp(as.matrix(t_c_matrix[]), center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(sample_bt)
autoplot(pca_object, data= t_c_matrix, colour = 'batch')
#
#remove batch effect-
adjusted <- ComBat_seq(as.matrix(counts_matrix.processed), batch=annotation$Batch, group=annotation$AD)
counts.norm <- adjusted/(colSums(adjusted)/1E6)

pca_object      <- prcomp(t(adjusted), center = TRUE, scale. = TRUE)
autoplot(pca_object, data= t_c_matrix, colour = 'batch')

counts_matrix.adj <- adjusted
#cpm_matrix.adj <-log2(counts_matrix.adj/(colSums(counts_matrix.adj)/1E6))

df_counts     <- as.data.frame(t(counts.norm))
df_counts$batch <- sample_bt
df_counts$AD <- sample_ID
counts2plot<- data.frame(counts = c(t(df_counts[,1:(ncol(df_counts)-2)])),batch = rep(df_counts[,ncol(df_counts)-1],nrow(counts.norm)),AD = rep(df_counts[,ncol(df_counts)],nrow(counts.norm)))#, AD = c(df_counts[,ncol(df_counts)]))
counts2plot$batch <- as.factor(counts2plot$batch)
counts2plot$counts <- as.numeric(counts2plot$counts)
counts2plot$AD <- as.factor(counts2plot$AD)
counts2plot$lcpm <- log10(counts2plot$counts)

p <- ggplot(counts2plot, aes(x=batch, y=lcpm)) +
  geom_boxplot() + theme_minimal() + scale_fill_wa_d(wacolors$volcano)#+ ylim(c(0,100))   #+
file_name <- paste(resultsPath,"/ROSMAP_cpm_boxplot.pdf",sep="")
ggsave(file_name, p, width = 12, height = 10, units = "cm",dpi = 400)

#get counts per million of transcripts
#counts.norm <- counts_matrix.processed/(colSums(adjusted)/1E6)
#for diff expr analysis⁄Keep those genes with at least 1 cpm in 50% of the NCI samples
AD_idxs <- which(annotation$AD=='AD')
NCIidxs <- which(annotation$AD=='No_AD')
nci <- counts.norm[,NCIidxs]
exprn <- rowSums((nci>=1))
toKeep <- which(exprn>=0.5*length(NCIidxs))
#Identify AD and non-AD smples
patient <- colnames(counts_matrix.processed)
annotation$patient <- patient

write.table(annotation,file = 'data/ROSMAP_annotation_processed.txt',row.names= FALSE,col.names=TRUE,sep= "\t")
tempMat <- adjusted[toKeep,]
tempMat <- counts_matrix.processed[toKeep,]
write.table(tempMat,file = 'data/ROSMAP_counts_processed_4DE.txt',row.names= FALSE,col.names=TRUE,sep= "\t")
genes <- rownames(adjusted)[toKeep]
write.table(genes,file = 'data/ROSMAP_genes_processed.txt',row.names= FALSE,col.names=TRUE,sep= "\t")
adjusted <- as.data.frame(adjusted)
adjusted$genes <- rownames(adjusted)
#get the pan-RNA sample
mat <- as.matrix(adjusted[,1:(ncol(adjusted)-1)])
#normalize to CPM and keep those genes with  >=1 CPM in at least one sample (permissive approach)
norm.counts <- mat/(colSums(mat)/1E6)
norm.counts.logical <- norm.counts>=1
norm.counts.logical <- rowSums(norm.counts.logical)
toKeep <- which(norm.counts.logical>=1)
adjusted <- adjusted[toKeep,]
write.table(adjusted,file = 'data/ROSMAP_counts_processed_all.txt',row.names= FALSE,col.names=TRUE,sep= "\t")
newDF <- data.frame(gene = adjusted$genes,counts = rowMeans(adjusted[,1:(ncol(adjusted)-1)]))
newDF$gene <- gsub("\\..*","",newDF$gene)
write.table(newDF,file = 'data/ROSMAP_meanCounts_processed_all.txt',row.names= FALSE,col.names=TRUE,sep= "\t")


#plot counts per batch after batch adjustment
# adjusted.f <- counts.norm[toKeep,]
# df_counts     <- as.data.frame(t(adjusted.f))
# df_counts$batch <- sample_bt
# df_counts$AD <- sample_ID
# counts2plot<- data.frame(counts = c(t(df_counts[,1:(ncol(df_counts)-2)])),batch = rep(df_counts[,ncol(df_counts)-1],nrow(adjusted.f)),AD = rep(df_counts[,ncol(df_counts)],nrow(adjusted.f)))#, AD = c(df_counts[,ncol(df_counts)]))
# counts2plot$batch <- as.factor(counts2plot$batch)
# counts2plot$counts <- as.numeric(counts2plot$counts)
# counts2plot$AD <- as.factor(counts2plot$AD)
# counts2plot$lcpm <- log10(counts2plot$counts)
# p <- ggplot(counts2plot, aes(x=batch, y=lcpm)) +
#   geom_boxplot() + theme_minimal() + scale_fill_wa_d(wacolors$volcano)#+ ylim(c(0,100))   #+
# file_name <- paste(resultsPath,"/ROSMAP_filt_cpm_boxplot.pdf",sep="")
# ggsave(file_name, p, width = 12, height = 10, units = "cm",dpi = 400)



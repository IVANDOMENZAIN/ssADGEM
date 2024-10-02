library("wacolors")
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("viridis")
library("plyr")
library('sva')
library("ggfortify")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/DE_analysis/",sep = "")
dir.create(resultsPath)
#load data
annotation <- read_delim(file = 'data/ROSMAP_annotation_processed.txt',delim = '\t', na='NA')
genes      <- read_delim(file = 'data/ROSMAP_genes_processed.txt',delim = '\t', na='NA')
counts     <- read_delim(file = 'data/ROSMAP_counts_processed_4DE.txt',delim = '\t', na='NA')
#danish_file<- read_delim(file = 'data/ROSMAP_counts_processed_4DE.txt',delim = '\t', na='NA')

# dgenes <- danish_file[,1]
# dcounts <- danish_file[,2:ncol(danish_file)]
# 
#  samples <- colnames(danish_file)
#  samples <- samples[2:length(samples)]
# #
# # #discard samples that are not in annotation
#  idxs <- match(samples,annotation$patient)
#  x <- which(is.na(idxs))
#  dcounts <- dcounts[,-x]
#  #now the other way around
#  idxs <- match(annotation$patient,samples)
#  x <- which(is.na(idxs))
#  annotation <- annotation[-x,]
# 
#  annotation %<>%
#    arrange(
#      match(patient, colnames(dcounts))
#    )
# counts <- dcounts

AD_idxs  <- which(annotation$AD=='AD')
NCIidxs  <- which(annotation$AD=='No_AD')

NCIdf <- counts[,NCIidxs]
sdevs <- rowSds(NCIdf)
meanV <- rowMeans(NCIdf)
coefV <- sdevs/meanV
toKeep <- which(coefV<1)
annotation$groups <- rep(0,nrow(annotation))
annotation$groups[AD_idxs] <- 1
annotation$groups[NCIidxs] <- 0
annotation$AD <- as.factor(annotation$AD)
levels(annotation$AD) <- c(1,0)
clusters <- c(1,2,3)

for (i in clusters){
cluster <- read_delim(file = paste('data/samples_GRN_cluster_',i,'.txt',sep=""),delim = ',', na='NA',skip = 0)
cluster <- colnames(cluster)
cluster <- gsub(' ','',cluster)
AD_idxs <- match(cluster,annotation$patient)
AD_idxs <- AD_idxs[!is.na(AD_idxs)]


#Perform DE analysis AD_i vs. lumped_NCI
#for (i in 1:length(AD_idxs)){
  newIdxs <- c(AD_idxs,NCIidxs)
  cmat <- as.data.frame(counts[toKeep,newIdxs])
  
  dataset <- DESeqDataSetFromMatrix(
    countData = cmat,
    colData = annotation[newIdxs,],
    design = ~AD
  )
  res  <- DESeq(dataset)
  res  <- results(res)
  DEdf <- data.frame(gene = genes$x[toKeep],meanVal  = res@listData$baseMean, log2FC = res@listData$log2FoldChange,pval = res@listData$pvalue,padj = res@listData$padj)
  DEdf$padj[is.na(DEdf$padj)] <- 1
  newDF <- DEdf[DEdf$padj<=0.01,]
  newDF <- newDF[order(-newDF$log2FC),]
  newDF$gene <- gsub("\\..*","",newDF$gene)
  DEdf$gene <- gsub("\\..*","",DEdf$gene)
  
  down <- which(DEdf$log2FC<0 & DEdf$padj<=0.01)
  up <- which(DEdf$log2FC>0 & DEdf$padj<=0.01)
  print(length(up))
    print(length(down))

  print(' ')
  down <- which(DEdf$log2FC<0 & DEdf$padj<=0.05)
  up <- which(DEdf$log2FC>0 & DEdf$padj<=0.05)
  print(length(down))
  print(length(up))
  write.table(DEdf,file = paste(resultsPath,'DE_genes_AD_GRN_cluster_',i,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")
  write.table(newDF,file = paste(resultsPath,'DE_genes_signif_AD_GRN_cluster_',i,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")
}
  
  
write.table(x,file = paste(resultsPath,'ivan_patients.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")

  
# #get a summary of DE results
# DE_genes <- data.frame(genes = genes,Dreg = rep(0,nrow(genes)),Ureg = rep(0,nrow(genes)))
# for (i in 1:length(which(annotation$AD=='AD'))){
#   DEdf <- read_delim(paste('results/DE_analysis/DE_genes_signif_AD_',i,'.txt',sep=""),delim= "\t")
#   up <- which(DEdf$log2FC>=1)
#   dn <- which(DEdf$log2FC<=-1)
#   upos <- match(DEdf$gene[up],genes$x)
#   dnos <- match(DEdf$gene[dn],genes$x)
#   DE_genes$Ureg[upos] <- DE_genes$Ureg[upos]+1
#   DE_genes$Dreg[dnos] <- DE_genes$Ureg[dnos]+1 
# }
# DE_genes$DEtotal <- (DE_genes$Ureg+DE_genes$Dreg)
# DE_genes <- DE_genes[order(-DE_genes$Ureg),]
# newDF <- DE_genes
# newDF$x <- gsub("\\..*","",newDF$x)
# write.table(newDF,file = 'results/ROSMAP_DE_genes_occurence.txt',row.names= FALSE,col.names=TRUE,sep= "\t")
# 
# genes$x <- gsub("\\..*","",genes$x)
# write.table(genes,file = 'data/ROSMAP_genes_4DE.txt',row.names= FALSE,col.names=TRUE,sep= "\t")

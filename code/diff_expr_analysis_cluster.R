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
annotation    <- read_delim(file = 'data/ROSMAP_annotation_processed.txt',delim = '\t', na='NA')
genes         <- read_delim(file = 'data/ROSMAP_genes_processed.txt',delim = '\t', na='NA')
counts         <- read_delim(file = 'data/ROSMAP_counts_processed_4DE.txt',delim = '\t', na='NA')

AD_idxs <- which(annotation$AD=='AD')
NCIidxs <- which(annotation$AD=='No_AD')

NCIdf <- counts[,NCIidxs]
sdevs <- rowSds(NCIdf)
meanV <- rowMeans(NCIdf)
coefV <- sdevs/meanV
toKeep <- which(coefV<1)
annotation$groups <- rep(0,nrow(annotation))
annotation$groups[AD_idxs] <- 2
annotation$groups[NCIidxs] <- 1
#Perform DE analysis AD_i vs. lumped_NCI
meanV <- meanV[toKeep]
for (i in 1:length(AD_idxs)){
  newIdxs <- c(AD_idxs[i],NCIidxs)
  cmat    <- as.data.frame(counts[toKeep,newIdxs])
  cmat    <- cmat/(colSums(cmat)/1E6)
  sample_DE <- data.frame(genes = genes[toKeep,],log2fc = rep(0,length(toKeep)),pvalue = rep(1,length(toKeep)))
  for (j in 1:length(toKeep)){
    vectorRef <- as.numeric(cmat[j,2:ncol(cmat)])
    observation <- as.numeric(cmat[j,1])
    metric <- ks.test(observation,observation,alternative = 'two.sided')
    sample_DE$pvalue[j] <- metric$p.value
    sample_DE$log2fc[j] <-log2(observation/mean(vectorRef))
  } 
  
  dataset <- DESeqDataSetFromMatrix(
    countData = cmat,
    colData = annotation[newIdxs,],
    design = ~groups
  )
  dds_p <- estimateSizeFactors(dataset)
  dds_p <- estimateDispersions(dds_p)
  res <- nbinomWaldTest(dds_p)
  res <- results(res)
  #dds <- DESeq(dataset)
  DEdf <- data.frame(gene = genes$x[toKeep],meanVal  = res@listData$baseMean, log2FC = res@listData$log2FoldChange,pval = res@listData$pvalue,padj = res@listData$padj)
  DEdf$pval[is.na(DEdf$padj)] <- 1
  newDF <- DEdf[abs(DEdf$log2FC) >=1 & DEdf$pval<=0.01,]
  newDF <- newDF[order(-newDF$log2FC),]
  write.table(DEdf,file = paste(resultsPath,'DE_genes_AD_',i,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")
  write.table(newDF,file = paste(resultsPath,'DE_genes_signif_AD_',i,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")
}
#get a summary of DE results
DE_genes <- data.frame(genes = genes,Dreg = rep(0,nrow(genes)),Ureg = rep(0,nrow(genes)))
for (i in 1:length(which(annotation$AD=='AD'))){
  DEdf <- read_delim(paste('results/DE_analysis/DE_genes_signif_AD_',i,'.txt',sep=""),delim= "\t")
  up <- which(DEdf$log2FC>=1)
  dn <- which(DEdf$log2FC<=-1)
  upos <- match(DEdf$gene[up],genes$x)
  dnos <- match(DEdf$gene[dn],genes$x)
  DE_genes$Ureg[upos] <- DE_genes$Ureg[upos]+1
  DE_genes$Dreg[dnos] <- DE_genes$Ureg[dnos]+1 
}
DE_genes$DEtotal <- (DE_genes$Ureg+DE_genes$Dreg)
DE_genes <- DE_genes[order(-DE_genes$Ureg),]
newDF <- DE_genes
newDF$x <- gsub("\\..*","",newDF$x)
write.table(newDF,file = 'results/ROSMAP_DE_genes_occurence.txt',row.names= FALSE,col.names=TRUE,sep= "\t")

genes$x <- gsub("\\..*","",genes$x)
write.table(genes,file = 'data/ROSMAP_genes_4DE.txt',row.names= FALSE,col.names=TRUE,sep= "\t")

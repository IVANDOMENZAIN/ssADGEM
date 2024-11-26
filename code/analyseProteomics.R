library("wacolors")
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("viridis")
library("plyr")
library('sva')
library("ggfortify")
library(tidyr)

#set wd and create the necessary ones for results
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/DE_prot_analysis/",sep = "")
dir.create(resultsPath)
source('code/translateuniprot.R')
#load clinical data 
annotation    <- read_delim(file = 'data/ROSMAP_annotation_processed.txt',delim = '\t', na='NA')
#load proteomics data and match sample IDs to those in the clinical data file (RXXXX IDs)
protData  <- read_delim(file = 'data/protemics_part1.txt',delim = '\t', na='NA')
colnames(protData) <- gsub('\\.','_',colnames(protData))
patients <- match(colnames(protData),annotation$proteomics_id)
patients <- patients[!is.na(patients)]
annotationProt <- annotation[patients,]

idxs <- match(annotation$proteomics_id[!is.na(annotation$proteomics_id)],colnames(protData))
idxs <- idxs[!is.na(idxs)]
proteins <- protData[,1]
protData <- protData[,idxs]
idxs <- match(colnames(protData),annotationProt$proteomics_id)
annotationProt <- annotationProt[idxs,]
patients <- annotationProt$patient
colnames(protData) <-patients
protData <- as.data.frame(protData)
rownames(protData) <- proteins$gene
x <- rownames(protData)
newDF <- as.data.frame(x)
newDF <- newDF %>%
  separate(x, into = c("prot_name", "uniprot"), sep = "\\|")
rownames(protData) <- newDF$uniprot
write.table(protData,file = 'data/proteomics_abundance_annotated.txt',row.names= TRUE,col.names=TRUE,sep= "\t")
write.table(annotationProt,file = 'data/annotation_proteomics.txt',row.names= FALSE,col.names=TRUE,sep= "\t")

#Discard proteins with NA values in more than 50% of samples
threshold  <- ncol(protData) * 0.5
protData.f <- protData[rowSums(is.na(protData)) < threshold, ]
#substitute remaining NA by zeros
protData.f[is.na(protData.f)] <- 0
#Discard proteins with 0 values in more than 50% of samples
protData.f <- protData.f[rowSums(protData.f == 0) <= threshold, ]
#plot PCA classified per batch
t_c_matrix       <- t(protData.f)
pca_object       <- prcomp(as.matrix(t_c_matrix), center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(annotationProt$proteomics_batch)
t_c_matrix$AD <- as.factor(annotationProt$AD)
autoplot(pca_object, data= t_c_matrix, colour = 'AD')
autoplot(pca_object, data= t_c_matrix, colour = 'batch')

#Outlier sample identified, remove 28_120411
x    <- pca_object$x
pos <- which(x[,1]==min(x[,1]))
protData.f <- protData.f[,-pos]
annotationProt <- annotationProt[-pos,]
#PCA again and Plot
t_c_matrix       <- t(protData.f)
pca_object       <- prcomp(as.matrix(t_c_matrix), center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(annotationProt$proteomics_batch)
t_c_matrix$AD <- as.factor(annotationProt$AD)
autoplot(pca_object, data= t_c_matrix, colour = 'batch')
autoplot(pca_object, data= t_c_matrix, colour = 'AD')

#Outlier sample identified, remove 04_120405
x    <- pca_object$x
pos <- which(x[,2]==max(x[,2]))
protData.f <- protData.f[,-pos]
annotationProt <- annotationProt[-pos,]
#PCA again and Plot
t_c_matrix       <- t(protData.f)
pca_object       <- prcomp(as.matrix(t_c_matrix), center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(annotationProt$proteomics_batch)
t_c_matrix$AD <- as.factor(annotationProt$AD)
autoplot(pca_object, data= t_c_matrix, colour = 'batch')
autoplot(pca_object, data= t_c_matrix, colour = 'AD')

#Outlier sample identified, remove "340_120501" "461_120514"
x    <- pca_object$x
pos <- which(x[,1]==min(x[,1]))
pos2 <- which(x[,2]==min(x[,2]))
toRmv <- c(pos,pos2)
protData.f <- protData.f[,-toRmv]
annotationProt <- annotationProt[-toRmv,]
#PCA again and Plot
t_c_matrix       <- t(protData.f)
pca_object       <- prcomp(as.matrix(t_c_matrix), center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$batch <- as.factor(annotationProt$proteomics_batch)
t_c_matrix$AD <- as.factor(annotationProt$AD)
t_c_matrix$ptfm <- as.factor(annotationProt$prot_platform)
t_c_matrix$tbtch <- as.factor(annotationProt$TMT_batch)

autoplot(pca_object, data= t_c_matrix, colour = 'batch')
autoplot(pca_object, data= t_c_matrix, colour = 'AD')
autoplot(pca_object, data= t_c_matrix, colour = 'ptfm')
autoplot(pca_object, data= t_c_matrix, colour = 'tbtch')
ADpos <- which(annotationProt$AD == 'AD')
NCIpos <- which(annotationProt$AD == 'No_AD')
annotationProt$groups <- '0'
annotationProt$groups[ADpos] <- '2'
annotationProt$groups[NCIpos] <- '1'
annotationProt$groups <- as.factor(annotationProt$groups)
#there are 55 samples remaining 20 AD and 35 NCI, diff. exp. by cluster
for (i in c(1,2,3)){
  cluster <- read_delim(file = paste('data/samples_GRN_cluster_',i,'.txt',sep=""),delim = ',', na='NA',skip = 0)
  cluster <- colnames(cluster)
  cluster <- gsub(' ','',cluster)
  pos <- match(cluster,annotationProt$patient)
  pos <- pos[!is.na(pos)]
  #subset data for DE analysis
  protdf <- protData.f[,c(pos,NCIpos)]
  metadf <- annotationProt[c(pos,NCIpos),]
  
  dataset <- DESeqDataSetFromMatrix(
    countData = trunc(protdf),
    colData = metadf,
    design = ~groups
  )
  res  <- DESeq(dataset)
  res  <- results(res)
  
  DEdf <- data.frame(protein = rownames(protdf),meanVal  = res@listData$baseMean, log2FC = res@listData$log2FoldChange,pval = res@listData$pvalue,padj = res@listData$padj)
  DEdf$padj[is.na(DEdf$padj)] <- 1
  newDF <- DEdf[DEdf$padj<=0.05,]
  newDF <- newDF[order(-newDF$log2FC),]
  newDF$protein <- gsub("\\..*","",newDF$protein)
  DEdf$protein <- gsub("\\..*","",DEdf$protein)
  
  DEdf <- DEdf %>%
    separate(protein, into = c("prot_name", "uniprot"), sep = "\\|")
  newDF <- newDF %>%
    separate(protein, into = c("prot_name", "uniprot"), sep = "\\|")
  
  down <- which(DEdf$log2FC<0 & DEdf$padj<=0.01)
  up <- which(DEdf$log2FC>0 & DEdf$padj<=0.01)
  print(length(up))
  print(length(down))
  
  print(' ')
  down <- which(DEdf$log2FC<0 & DEdf$padj<=0.05)
  up <- which(DEdf$log2FC>0 & DEdf$padj<=0.05)
  print(length(up))
  print(length(down))
  DEdf <- translateuniprot(DEdf,uniprot)
  write.table(DEdf,file = paste(resultsPath,'DE_proteins_AD_GRN_cluster_',i,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")
  newDF <- translateuniprot(newDF,uniprot)
  write.table(newDF,file = paste(resultsPath,'DE_proteins_signif_AD_GRN_cluster_',i,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")

  
}
  
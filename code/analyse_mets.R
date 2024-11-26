#Analyse metabolomics
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
resultsPath <- paste(resultsPath,"/metabolomics/",sep = "")
dir.create(resultsPath)
#load data
annotation <- read_delim(file = 'data/ROSMAP_annotation_processed_clusters.txt',delim = '\t', na='NA')
data      <- read_delim(file = 'data/ROSMAP_metabolomics_data.txt',delim = '\t', na='NA')
metIDs <- read_delim(file = 'data/data_dictionary.txt',delim = '\t', na='NA')

#REduce dataset to those samples present in annotation (from RNAseq studies)
idxs <- match(data$individualID,annotation$individualID)
toKeep <- which(!is.na(idxs))
sampleIDs <- data$individualID[toKeep]
idxs <- match(annotation$individualID,sampleIDs)
toKeep <- which(!is.na(idxs))
annotation <- annotation[toKeep,]
annotation$cluster[which(is.na(annotation$cluster))]<- 0
data.df <- as.data.frame(data[toKeep,2:ncol(data)])
data.na <-  is.na(data.df)
data.na <- t(data.na)
data.na <- rowSums(data.na)
#Keep just those columns wiht NA values in less than 50% of the samples
toRmv <- which(data.na>0.5*nrow(data.df))
data.df <- data.df[,-toRmv]
newDF <- data.df
chemicalNames <- colnames(newDF)
idxs <- match(chemicalNames,metIDs$CHEM_ID)
#substitute NA values by lowest read in the sample for each metabolite
for (i in 1:ncol(data.df)){
 
  nonNa <- data.df[!is.na(data.df[,i]),i]
  minval <- min(nonNa)
  newDF[which(is.na(data.df[,i])),i] <- 0.5*minval
}

idxs <- match(annotation$individualID,sampleIDs)
newDF <- newDF[idxs,]
sampleIDs <- sampleIDs[idxs]
#ASsess global differences
t_c_matrix       <- newDF
#t_c_matrix$AD <- 
pca_object       <- prcomp(t_c_matrix, center = TRUE, scale. = TRUE)
t_c_matrix       <- as.data.frame(t_c_matrix)
t_c_matrix$cluster <- as.factor(annotation$cluster)
autoplot(pca_object, data= t_c_matrix, colour = 'cluster')
#assess significance met by met
common <- c()
for(j in 1:3){
  cluster <- read_delim(file = paste('data/samples_GRN_cluster_',j,'.txt',sep=""),delim = ',', na='NA',skip = 0)
  cluster <- colnames(cluster)
  cluster <- gsub(' ','',cluster)
  AD_idxs <- match(cluster,annotation$patient)
  AD_idxs <- AD_idxs[!is.na(AD_idxs)] 
  
  significant <- c()
  results<- data_frame(foldchange= c(),pvalue=c())
  fchanges <- c()
  pvalues <- c()
  for (i in 1:ncol(newDF)){
    ADdist  <- newDF[AD_idxs,i]
    meanAD <- mean(ADdist)
    NCIdist <- newDF[which(annotation$AD=='No_AD'),i]
    meanNCI <- mean(NCIdist)
    result <- wilcox.test(ADdist, NCIdist, exact = TRUE)
    fchange <- meanAD/meanNCI
    fchanges <- c(fchanges,fchange)
    pvalues  <- c(pvalues,result$p.value)
    # 
    #   if (result$p.value<0.05 & fchange<1){
    #     chemical <- colnames(newDF)[i]
    #     idx <- which(metIDs$CHEM_ID==chemical)
    #     chemName <- metIDs$SHORT_NAME[idx]
    #     significant <- c(significant,chemName)
    #     print(chemName)
    #   }
  }
  print(' ')
  adjpvalues <- p.adjust(pvalues, method = "BH")
  sigMets <- data.frame(idxs =colnames(newDF),fchanges = fchanges, pvalues=pvalues,adjpvalues=adjpvalues)
  idxs1 <- which(pvalues<=0.05)
  if (length(idxs>0)){
    chemicals <- colnames(newDF)[idxs1]
    idxs <- match(chemicals,metIDs$CHEM_ID)
    chemNames <- metIDs$SHORT_NAME[idxs]
    print(chemNames)
  }
  resultsDF <- data.frame(chemNames,metIDs$SUPER_PATHWAY[idxs],metIDs$SUB_PATHWAY[idxs],metIDs$KEGG[idxs],fchanges[idxs1],pvalues[idxs1],adjpvalues[idxs1])
  write.table(resultsDF,file = paste(resultsPath,'diff_mets_cluster_',j,'.txt',sep=""),row.names= FALSE,col.names=TRUE,sep= "\t")
  
}




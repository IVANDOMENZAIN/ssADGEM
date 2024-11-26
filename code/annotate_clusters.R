library("tidyverse") # Tibble dataframes
library("plyr")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()

#load data
annotation <- read_delim(file = 'data/ROSMAP_annotation_processed.txt',delim = '\t', na='NA')
genes      <- read_delim(file = 'data/ROSMAP_genes_processed.txt',delim = '\t', na='NA')
counts     <- read_delim(file = 'data/ROSMAP_counts_processed_4DE.txt',delim = '\t', na='NA')
metData    <- read_delim(file = 'data/ROSMAP_metabolomics_data.txt',delim = '\t', na='NA')
metIDs <- metData$individualID
overlap <- intersect(metIDs,annotation$individualID)


NCI <- annotation$individualID[which(annotation$AD == 'No_AD')]
overlap <- intersect(metIDs,NCI)
#danish_file<- read_delim(file = 'data/ROSMAP_counts_processed_4DE.txt',delim = '\t', na='NA')
annotation$cluster <- rep(NA,length(annotation$individualID))
# dgenes <- danish_file[,1]
for (i in c(1,2,3)){
  cluster <- read_delim(file = paste('data/samples_GRN_cluster_',i,'.txt',sep=''),delim = ', ', na='NA')
  cluster <- colnames(cluster)
  idxs    <- match(cluster,annotation$specimenID)
  idxs    <- idxs[!is.na(idxs)]
  annotation$cluster[idxs] <- i
  
  RNA <- annotation$individualID[idxs]
  overlap <- intersect(metIDs,RNA)
  print(length(overlap))
}
write_delim(annotation, file = 'data/ROSMAP_annotation_processed_clusters.txt',delim = '\t', na='NA')

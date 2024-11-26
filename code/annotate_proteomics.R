annotate_proteomics <- function(annotation){

#load data
#annotation    <- read_delim(file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')

protData  <- read_delim(file = 'data/protemics_part1.txt',delim = '\t', na='NA')
metadata  <- read_delim(file = 'data/ROSMAP_assay_proteomics_TMTquantitation_metadata.csv',delim = ',', na='NA')
metadata2 <- read_delim(file = 'data/ROSMAP_biospecimen_metadata.csv',delim = ',', na='NA')
clinical  <- read_delim(file = 'data/ROSMAP_clinical.csv',delim = ',', na='NA')
#create a proteomics ID field in annotation df
annotation$proteomics_id <- rep(NA,length(annotation$individualID))
annotation$proteomics_batch <- rep(NA,length(annotation$individualID))
annotation$TMT_batch <- rep(NA,length(annotation$individualID))
annotation$prot_platform <- rep(NA,length(annotation$individualID))
patients <- paste('.',annotation$individualID,sep='')
i <- 1
for (patient in patients){
  idxs <- grep(patient,metadata$specimenID)
  if (length(idxs)>0){
    annotation$proteomics_batch[i] <-  metadata$batch[idxs[1]]
    annotation$TMT_batch  <-  metadata$TMTdataSubmissionBatch[idxs[1]]
    annotation$proteomics_id[i] <- metadata$batchChannel[idxs[1]]
    annotation$prot_platform[i] <- metadata$platform[idxs[1]]
  }
  i <- i +1 
}
annotation$proteomics_id <- gsub('\\.','_',annotation$proteomics_id)
#update clinical data file
#write_delim(annotation, file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')
return(annotation)
}
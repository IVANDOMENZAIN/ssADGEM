# Installing required R packages ===============================================
if (!requireNamespace("rstudioapi", quietly = TRUE)){
  install.packages("rstudioapi")}
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")}
if (!requireNamespace("magrittr", quietly = TRUE)){
  install.packages("magrittr")}
if (!requireNamespace("DESeq2", quietly = TRUE)){
  install.packages("DESeq2")}
if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
# Changing working directory ===================================================
# Setting the working directory to the directory containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Needed libraries
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
#data 
source('load_data.R')
output               <- load_data()
count_matrix         <-output[[1]]
rnaseq_metadata      <- output[[2]]
clinical_metadata    <- output[[3]]
biospecimen_metadata <- output[[4]]

sample_ids   <- colnames(count_matrix)
gene_ids     <- rownames(count_matrix)[5:nrow(count_matrix)]
count_matrix <- count_matrix[5:nrow(count_matrix),]
mode(count_matrix) <- "integer"
# Merge annotation data into one dataframe by specimenID and individualID
annotation_df <- merge(rnaseq_metadata, biospecimen_metadata,
                       by = "specimenID",
                       all = TRUE)

annotation_df %<>% merge(clinical_metadata,
                         by = "individualID",
                         all = TRUE)

## annotation/count_matrix has some entries that should be removed

### Duplicate entry
if (all(count_matrix[, "150_120419"] == count_matrix[, "150_120419_0_merged"])) {
  idx <- which(dimnames(count_matrix)[[2]] == "150_120419_0_merged")
  count_matrix <- count_matrix[, -idx]
}

annotation_df %<>%
  filter(specimenID %in% colnames(count_matrix)) %>%
  filter(assay == "rnaSeq") %>%
  filter(!is.na(cogdx)) %>%
  filter(!is.na(braaksc)) %>%
  filter(!is.na(ceradsc)) %>%
  # RINcontinous is missing for ~ 100 entries, but RIN is only missing one 
  mutate(RINcontinuous = coalesce(RINcontinuous, RIN)) %>%
  filter(!is.na(RINcontinuous)) %>%
  filter(!is.na(pmi)) %>%
  filter(!is.na(age_death)) #%>%
  #filter((exclude == FALSE) | is.na(exclude))

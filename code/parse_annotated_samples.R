# # Installing required R packages ===============================================
# if (!requireNamespace("rstudioapi", quietly = TRUE)){
#   install.packages("rstudioapi")}
# if (!requireNamespace("tidyverse", quietly = TRUE)){
#   install.packages("tidyverse")}
# if (!requireNamespace("magrittr", quietly = TRUE)){
#   install.packages("magrittr")}
# if (!requireNamespace("DESeq2", quietly = TRUE)){
#   install.packages("DESeq2")}
# if (!require("DESeq2", quietly = TRUE))
#   BiocManager::install("DESeq2")
# # Changing working directory ===================================================
# # Setting the working directory to the directory containing this repository
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd('..')
# 
# 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()

# Needed libraries
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
#data 
source('code/load_data.R')
output               <- load_data()
count_matrix         <-output[[1]]
rnaseq_metadata      <- output[[2]]
clinical_metadata    <- output[[3]]
biospecimen_metadata <- output[[4]]
sample_ids   <- colnames(count_matrix)
gene_ids     <- rownames(count_matrix)[5:nrow(count_matrix)]
#count_matrix <- count_matrix[5:nrow(count_matrix),]
mode(count_matrix) <- "integer"
# Merge annotation data into one dataframe by specimenID and individualID
annotation_df <- merge(rnaseq_metadata, biospecimen_metadata,
                       by = "specimenID",
                       all = TRUE)

annotation_df %<>% merge(clinical_metadata,
                         by = "individualID",
                         all = TRUE)
## annotation/count_matrix has some entries that should be removed

### REmove Duplicate entry
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
  #mutate(RINcontinuous = coalesce(RINcontinuous, RIN)) %>%
  filter(!is.na(RIN)) %>%
  filter(!is.na(pmi)) %>%
  filter(!is.na(age_death)) %>%
  filter((exclude == FALSE) | is.na(exclude))
## Some of the annotation columns are mostly NA values, remove such columns
annotation_df %<>%
  select_if(function(col)
    n_distinct(col) > 1 &
      sum(!is.na(col)) > .9 * length(col) #more than 90% are non-na
  )
## Some rows/columns are now redundant
annotation_df %<>% distinct
annotation_df <- annotation_df[, !duplicated(annotation_df %>% t)]
## Sample 492_120515 belongs to 2 batches,
## assume the latter one is correct
annotation_df %<>% filter(specimenID != "492_120515_0" | is.na(specimenID))
## Remove the counts vectors lacking metadata (remove from the counts matrix)
count_matrix <- count_matrix[, annotation_df$specimenID]
## counts matrix contains global info about each sample (#mapped, #unmapped, etc)
## transfer this info to metadata
idx <- count_matrix %>%
  dimnames %>%
  .[[1]] %>%
  grep("^N_*", .)

rownames(annotation_df) <- annotation_df$specimenID
#
annotation_df %<>%
  merge(count_matrix[idx,] %>% t %>% as.data.frame,
        by = 'row.names') %>%
  mutate(Row.names = NULL)

count_matrix <- count_matrix[-idx,]
## Sort annotations to match order in counts
annotation_df %<>%
  arrange(
    match(specimenID, colnames(count_matrix))
  )
rownames(annotation_df) <- annotation_df$specimenID # merge removes the rownames
# Make sure that counts and annotations have the same order
is_identical <- rownames(annotation_df) == colnames(count_matrix)
if (!all(is_identical)) {
  stop("There are unannotated samples in the counts", call. = FALSE)
}
# Change numeric codes to factors, add explanatory
# level names to some, rename columns, and remove unneeded columns
## Use only rownames (specimenID) as keys

# annotation_df %<>%
#   mutate(
#     specimenID = NULL,
#     projid = NULL,
#     individualID = NULL,
#     ID = NULL,
#     Sampleid = NULL
#   )

## Combine information from libraryBatch and Batch
annotation_df$libraryBatch %<>% as.factor
#rename some columns
names(annotation_df)[names(annotation_df) == "libraryBatch"] <- "Batch"
names(annotation_df)[names(annotation_df) == "spanish"] <- "latinx"
annotation_df$Study %<>% as.factor
annotation_df$msex %<>% as.factor
levels(annotation_df$msex) <- c("Female", "Male")
levels(annotation_df$Study) <- c("MAP", "ROS")
annotation_df$race %<>% as.factor
annotation_df$latinx %<>% as.factor
annotation_df$apoe_genotype %<>% as.factor
annotation_df$braaksc %<>% as.factor
annotation_df$cogdx %<>% as.factor
levels(annotation_df$cogdx) <- c("NCI", "MCI", "MCI+", "AD", "AD+", "Other")
annotation_df$dcfdx_lv %<>% as.factor
levels(annotation_df$dcfdx_lv) <- c("NCI", "MCI", "MCI+", "AD", "AD+", "Other")
## Add a binary ceradsc
annotation_df$ceradsc_binary <- as.numeric(annotation_df$ceradsc)
annotation_df$ceradsc %<>% as.factor
levels(annotation_df$ceradsc) <- c("Definite", "Probable", "Possible", "No_AD")
annotation_df$ceradsc_binary[annotation_df$ceradsc_binary== 1] <- "AD"
annotation_df$ceradsc_binary[annotation_df$ceradsc_binary== 2] <- "other"
annotation_df$ceradsc_binary[annotation_df$ceradsc_binary== 3] <- "other"
annotation_df$ceradsc_binary[annotation_df$ceradsc_binary== 4] <- "No_AD"
#annotation_df$ceradsc_binary %<>% as.factor

## Change age from string to numeric
annotation_df$age_at_visit_max %<>%
  sub("\\+", "", .) %>% # 90+ is treated as 90
  as.numeric
annotation_df$age_death %<>%
  sub("\\+", "", .) %>% # 90+ is treated as 90
  as.numeric
annotation_df$braaksc %<>% as.numeric
#Establish a definite and stringent AD diagnosis
ADpos <- which((annotation_df$ceradsc_binary == 'AD'&annotation_df$cogdx == 'AD'))
NoADpos <- which((annotation_df$ceradsc == 'No_AD' & annotation_df$cogdx == 'NCI'))
annotation_df$AD <- rep('other',nrow(annotation_df))
annotation_df$AD[ADpos] <- 'AD'
annotation_df$AD[NoADpos] <- 'No_AD'

gene_ids <- rownames(count_matrix)
matC <- as.data.frame(count_matrix)
rownames(matC) <- gene_ids
gene_ids <- as.data.frame(gene_ids)
write_delim(gene_ids, file = 'data/ROSMAP_annotated_samples_geneIDs.txt',delim = '\t', na='NA')
write_delim(matC, file = 'data/ROSMAP_annotated_samples_counts.txt',delim = '\t', na='NA')
write_delim(annotation_df, file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')
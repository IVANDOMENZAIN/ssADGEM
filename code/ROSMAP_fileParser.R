#!/usr/bin/env Rscript

# Parse arguments from terminal
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Format: ./<script> <data-folder> <target>", call. = FALSE)
  }
  directory_path <- args[1]
  target_path <- args[2]
} else {
  # Bad code, remember to remove
  directory_path <- "/castor/project/home/gasto/sens2023522/synapseFiles"
}


# Needed libraries
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")

# Read all files
working_directory <- getwd()
setwd(directory_path)

count_matrix <- "ROSMAP_all_counts_matrix.txt" %>%
  read_delim(delim = "\t") %>%
  column_to_rownames("feature") %>%
  as.matrix()

rnaseq_metadata <- read_csv("ROSMAP_assay_rnaSeq_metadata.csv")
clinical_metadata <- read_csv("ROSMAP_clinical.csv")
biospecimen_metadata <- read_csv("ROSMAP_biospecimen_metadata.csv")

setwd(working_directory)

# Merge annotation data into one dataframe by specimenID and individualID
annotation_df <- merge(rnaseq_metadata, biospecimen_metadata,
                       by = "specimenID",
                       all = TRUE
                     )
annotation_df %<>% merge(clinical_metadata,
                         by = "individualID",
                         all = TRUE
                       )

# Modify to align row and column names

## Format in file is 'X<specimenID>', we substitute 'X' to ''
dimnames(count_matrix)[[2]] %<>% sub("X", "", .)

## count_matrix has some that should be removed
if (all(count_matrix[, "150_120419"] == count_matrix[, "150_120419_0_merged"])) {
  
  idx <- which(dimnames(count_matrix)[[2]] == "150_120419_0_merged")
  count_matrix <- count_matrix[, -idx]
}


for (col_id in count_matrix %>% colnames()) {
  
  if (col_id %in% biospecimen_metadata$specimenID) {
    idx <- biospecimen_metadata %>%
      nrow() %>% 1:. %>%
      .[biospecimen_metadata$specimenID == col_id]
    
    is_excluded <- idx %>%
      biospecimen_metadata$exclude[.] %>%
      any(na.rm = TRUE)
    
    if (is_excluded) {
      which() STARTHERE
      count_matrix <- count_matrix[, -col_id]
    }
  }
}

## We can pick out rows of interest and discard the rest,
## additionally we sort them to the same order
id_of_interest <- count_matrix %>% colnames()

### By rows
annotation_df %<>%
  filter(
    specimenID %in% id_of_interest, # Corresponds to counts
    assay == "rnaSeq", # Is RNA-Seq
    exclude %>% is.na() # Not marked as excluded
  )

### By columns
annotation_df %<>%
  select_if(
    ~ n_distinct(.) > 1 # Has different values
  )

### Sort to match order in counts
annotation_df %<>%
  arrange(
    match(specimenID, id_of_interest)
  )

### Some are missing!!!
# Two samples were swapped in biospec and are excluded
# There is a redo and merge???
missing <- biospecimen_metadata %>%
  filter(
    specimenID %in% id_of_interest,
    !(specimenID %in% annotation_df$specimenID)
  )

# Construct the DESeq data set (dds), no design in this function
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = annotation_df,
  design = ~ "ADD DESIGN FOR AD OR NOT"
  )

# Save the dds
save(dds, file = target_path)
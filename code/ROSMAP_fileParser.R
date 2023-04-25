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
  # For using interactively, change paths to appropriate
  directory_path <- "/castor/project/home/gasto/sens2023522/synapseFiles"
  target_path <- "/castor/project/home/gasto/ssADGEM/data/ROSMAP_dds.rds"
}


# Needed libraries
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")

# Read all files
working_directory <- getwd()
directory_path %>% setwd

count_matrix <- "ROSMAP_all_counts_matrix.txt" %>%
  read_delim(delim = "\t") %>%
  column_to_rownames("feature") %>%
  as.matrix()

rnaseq_metadata <- read_csv("ROSMAP_assay_rnaSeq_metadata.csv")
clinical_metadata <- read_csv("ROSMAP_clinical.csv")
biospecimen_metadata <- read_csv("ROSMAP_biospecimen_metadata.csv")

working_directory %>% setwd

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

  if (col_id %in% annotation_df$specimenID) {

    idx <- which(
      annotation_df$specimenID == col_id &
      annotation_df$assay == "rnaSeq"
      )

    is_excluded <- idx %>%
      annotation_df$exclude[.] %>%
      any(na.rm = TRUE)

    lacks_cogdx <- idx %>%
      annotation_df$cogdx[.] %>%
      is.na

    if (is_excluded || lacks_cogdx) {
      col_idx <- which(
        dimnames(count_matrix)[[2]] == col_id
        )
      count_matrix <- count_matrix[, -col_idx]
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

# Make sure there are no missing values
missing <- biospecimen_metadata %>%
  filter(
    specimenID %in% id_of_interest,
    !(specimenID %in% annotation_df$specimenID)
  )
if (nrow(missing) > 0) {
  stop("There are unannotated samples in the counts", call. = FALSE)
}


# Construct the DESeq data set (dds),
# design is irrelevant as it will be changed later
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = annotation_df,
  design = ~ 1
  )

# Save the dds
target_path %>%
  dirname %>%
  setwd

target_path %>%
  basename %>%
  save(dds, file = ., compress = TRUE)

working_directory %>% setwd
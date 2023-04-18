#!/usr/bin/env Rscript

# Parse arguments from terminal
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Format: ./<script> <data-folder> <target>", call. = FALSE)
}
directory_path <- args[1]
target_path <- args[2]

# Needed libraries
library("DESeq2")

# Read all files
working_directory <- getwd()
setwd(directory_path)

count_matrix <- as.matrix(read.csv(
  "ROSMAP_all_counts_matrix.txt",
  sep = "\t",
  row.names = "feature"
  ))
rnaseq_metadata <- read.csv("ROSMAP_assay_rnaSeq_metadata.csv")
clinical_metadata <- read.csv("ROSMAP_clinical.csv")
biospecimen_metadata <- read.csv("ROSMAP_biospecimen_metadata.csv")

setwd(working_directory)

# Modify to align row and column names

## Format in file is 'X<specimenID>', we substitute 'X' to ''
dimnames(count_matrix)[[2]] <- sub("X", "", dimnames(count_matrix)[[2]])

## Function to remove na columns and rows
clean_df <- function(df) {
  # Remove na cols
  good_cols <- c()
  for (col_idx in seq_along(length(df))) {
    if (any(complete.cases(df[[col_idx]])))
      good_cols <- append(good_cols, col_idx)
  }
  df <- df[, good_cols]
  # Remove na rows
  good_rows <- c()
  for (row_idx in seq_along(nrow(df))) {
    if (any(complete.cases(df[row_idx, ])))
      good_rows <- append(good_rows, row_idx)
  }
  df <- df[good_rows, ]
  # Remove cols with no info
  good_cols <- c()
  for (col_idx in seq_along(length(df))) {
    if (!(length(unique(df[[col_idx]])) == 1))
      good_cols <- append(good_cols, col_idx)
  }
  df <- df[, good_cols]
  return(df)
}


## We can pick out rows of interest and discard the rest,
## additionally we sort them to the same order
id_of_interest <- colnames(count_matrix)

idx <- match(id_of_interest, rnaseq_metadata$specimenID)
rnaseq_metadata <- rnaseq_metadata[idx, ]
rnaseq_metadata <- cleanDF(rnaseq_metadata)

idx <- match(id_of_interest, biospecimen_metadata$specimenID)
biospecimen_metadata <- biospecimen_metadata[idx, ]
biospecimen_metadata <- cleanDF(biospecimen_metadata)

### Clinical metadata has no specimenID
individual_id_of_interest <- biospecimen_metadata$individualID
idx <- match(individual_id_of_interest, clinical_metadata$individualID)
clinical_metadata <- clinical_metadata[idx, ]

# Merge them into one dataframe by specimenID
annotation_df <- merge(rnaseq_metadata, biospecimen_metadata, by = "specimenID")
annotation_df <- merge(annotation_df, clinical_metadata, by = "individualID")

#
# Add na rows to counts with no annotation
#

# Construct the DESeq data set (dds), no design in this function
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = annotation_df,
  design = ~ 1
  )

# Save the dds
save(dds, file = target_path)
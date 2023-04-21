#!/usr/bin/env Rscript

# Parse arguments from terminal
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Format: ./<script> <data-folder> <target>", call. = FALSE)
}
directory_path <- args[1]
target_path <- args[2]

# Needed libraries
library("tidyverse")
library("DESeq2")

# Read all files
working_directory <- getwd()
setwd(directory_path)

count_matrix <- as.matrix(read_csv(
  "ROSMAP_all_counts_matrix.txt",
  sep = "\t",
  row.names = "feature"
  ))
rnaseq_metadata <- read_csv("ROSMAP_assay_rnaSeq_metadata.csv")
clinical_metadata <- read_csv("ROSMAP_clinical.csv")
biospecimen_metadata <- read_csv("ROSMAP_biospecimen_metadata.csv")

setwd(working_directory)

# Merge annotation data into one dataframe by specimenID and individualID
annotation_df <- merge(rnaseq_metadata, biospecimen_metadata, by = "specimenID")
annotation_df %<>% merge(clinical_metadata, by = "individualID")

# Modify to align row and column names

## Format in file is 'X<specimenID>', we substitute 'X' to ''
dimnames(count_matrix)[[2]] %<>% sub("X", "", .)

## We can pick out rows of interest and discard the rest,
## additionally we sort them to the same order
id_of_interest <- count_matrix %>% colnames()

annotation_df %>% filter(any(id_of_interest == specimenID))
idx <- match(id_of_interest, rnaseq_metadata$specimenID)
rnaseq_metadata <- rnaseq_metadata[idx, ]
rnaseq_metadata <- cleanDF(rnaseq_metadata)


## Function to remove na columns and rows
clean_df <- function(df) {

  # Remove na cols
  informative_col <- c()
  for (col_idx in seq_along(length(df))) {

    if (any(complete.cases(df[[col_idx]])))
      informative_col <- append(informative_col, col_idx)
  }
  df <- df[, informative_col]
  
  # Remove na rows
  informative_rows <- c()
  for (row_idx in seq_along(nrow(df))) {
  
    if (any(complete.cases(df[row_idx, ])))
      informative_rows <- append(informative_rows, row_idx)
  }
  df <- df[informative_rows, ]
  
  # Remove cols with only one value
  informative_col <- c()
  for (col_idx in seq_along(length(df))) {
  
    if (!(length(unique(df[[col_idx]])) == 1))
      informative_col <- append(informative_col, col_idx)
  }
  df <- df[, informative_col]
  
  return(df)
}
# pseudo
# get idx of all specimenID not in annotation
# missing_specID = func(count_matrix, annotation_df)
# 
# for( specID in missing ) {
#   add a row to data frame with only specimenID 
# }
#
#
#
# Add na rows to counts with no annotation
# Maybe for loop through, insert row at each discongru?
# Possibly another library has a good function for this
#

# Construct the DESeq data set (dds), no design in this function
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = annotation_df,
  design = ~ "ADD DESIGN FOR AD OR NOT"
  )

# Save the dds
save(dds, file = target_path)
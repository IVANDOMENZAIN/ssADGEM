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

count_matrix <- as.matrix(read.csv("ROSMAP_all_counts_matrix.txt"))
rnaseq_metadata <- read.csv("ROSMAP_assay_RNAseq_metadata.csv")
clinical_metadata <- read.csv("ROSMAP_clinical.csv")
biospecimen_metadata <- read.csv("ROSMAP_biospecimen_metadata.csv")

setwd(working_directory)

# Modify to align row and column names

## Format in file is 'X<specimenID>', we substitute 'X' to ''
colnames(count_matrix) <- sub("X", "", colnames(count_matrix))

## Clinical metadata has no specimenID, we check corresponding
## specimenID from biospecimen metadata
idx <- match(
    clinical_metadata$individualID,
    biospecimen_metadata$individualID
    )
clinical_metadata$specimenID <- biospecimen_metadata$specimenID[idx]

## We can pick out rows of interest and discard the rest,
## additionally we sort them to the same order
id_of_interest <- colnames(count_matrix)

tmp <- match(id_of_interest, rnaseq_metadata$specimenID)
rnaseq_metadata <- rnaseq_metadata[, colnames(count_matrix)]

tmp <- match(id_of_interest, clinical_metadata$specimenID)
clinical_metadata <- clinical_metadata[, colnames(count_matrix)]

tmp <- match(id_of_interest, biospecimen_metadata$specimenID)
biospecimen_metadata <- biospecimen_metadata[, colnames(count_matrix)]

# Merge them into one dataframe by specimenID
annotation_df <- merge(rnaseq_metadata, clinical_metadata, by = "specimenID")
annotation_df <- merge(annotation_df, biospecimen_metadata, by = "specimenID")

# Construct the DESeq data set (dds), no design in this function
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = annotation_df,
  design = ~ 1
  )

# Save the dds
save(dds, file = target_path)

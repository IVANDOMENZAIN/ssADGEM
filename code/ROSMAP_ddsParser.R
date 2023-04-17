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

## We can discard rows that have no column in counts
###HOW DO I DO THIS??? INVESTIGATE
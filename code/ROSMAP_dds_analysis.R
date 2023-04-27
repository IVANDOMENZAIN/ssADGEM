#!/usr/bin/env Rscript

# Parse arguments from terminal
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Format: ./<script> <dds-file> <target>", call. = FALSE)
  }
  dds_path <- args[1]
  target_path <- args[2]
} else {
  # For using interactively, change paths to appropriate
  dds_path <- "/castor/project/home/gasto/ssADGEM/data/ROSMAP_dds.rds"
  target_path <- NA
}

# Needed libraries
library("magrittr") # Piping
library("DESeq2")

# DE analysis
dds %<>% DESeq
dds %>% results(contrast = c(cogdx, "4", "1"))
dds %>% resultsNames
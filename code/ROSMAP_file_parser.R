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
  project_path <- "/castor/project/home/gasto/ssADGEM/"
  directory_path <- paste(project_path, "external/synapse_dir/", sep = "")
  target_path <- paste(project_path, "data/ROSMAP_dds.rds", sep = "")
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

tech_metadata <- "AMP_AD_ROSMAP_Broad-Rush_RNA_Seq_RIN.txt" %>%
  read.table(sep='\t',header=T)

picard_metadata <- "ROSMAP_all_metrics_matrix.txt" %>%
  read_delim(delim = "\t")
picard_metadata$specimenID <- picard_metadata$sample

## Unsure here
picard_important_cols <- c("specimenID",
                           "AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED",
                           "RnaSeqMetrics__PCT_CODING_BASES",
                           "RnaSeqMetrics__PCT_INTERGENIC_BASES",
                           "RnaSeqMetrics__PCT_INTRONIC_BASES", 
                           "RnaSeqMetrics__PCT_RIBOSOMAL_BASES")
picard_metadata %<>% .[, picard_important_cols]

working_directory %>% setwd

# Change count data from double to integer (smaller object)
mode(count_matrix) <- "integer"

# Merge annotation data into one dataframe by specimenID and individualID
annotation_df <- merge(rnaseq_metadata, biospecimen_metadata,
                       by = "specimenID",
                       all = TRUE)

annotation_df %<>% merge(clinical_metadata,
                         by = "individualID",
                         all = TRUE)

annotation_df %<>% merge(tech_metadata,
                         by = "projid",
                         all = TRUE)

annotation_df %<>% merge(picard_metadata,
                         by = "specimenID",
                         all = TRUE)

# Modify count and annotations to remove bad/unnecessary data,
# also align row and column names

## Format in file is 'X<specimenID>', we substitute 'X' to ''
dimnames(count_matrix)[[2]] %<>% sub("X", "", .)

## annotation/count_matrix has some entries that should be removed

### Duplicate entry
if (all(count_matrix[, "150_120419"] == count_matrix[, "150_120419_0_merged"])) {
  idx <- which(dimnames(count_matrix)[[2]] == "150_120419_0_merged")
  count_matrix <- count_matrix[, -idx]
}

### Some samples are lacking critical metadata
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
  filter(!is.na(age_death)) %>%
  filter((exclude == FALSE) | is.na(exclude))

## Some of the annotation columns have no variation
annotation_df %<>%
  select_if(function(col)
    n_distinct(col) > 1 &
    count(!is.na(col)) > .9 * length(col) #more than 90% are non-na
  )

## Some rows/columns are now redundant
annotation_df %<>% distinct
annotation_df <- annotation_df[, !duplicated(annotation_df %>% t)]

## Sample 492_120515 belongs to 2 batches,
## assume the latter one is correct
annotation_df %<>% filter(ID != "492_120515_0" | is.na(ID))

## Remove the counts lacking metadata
count_matrix <- count_matrix[, annotation_df$specimenID]

### Change rows for number of unmapped, etc. into metadata
idx <- count_matrix %>%
  dimnames %>%
  .[[1]] %>%
  grep("^N_*", .)
rownames(annotation_df) <- annotation_df$specimenID
annotation_df %<>%
  merge(count_matrix[idx,] %>% t %>% as.data.frame,
                                by = 'row.names') %>%
  mutate(Row.names = NULL)
rownames(annotation_df) <- annotation_df$specimenID # merge removes the rownames
count_matrix <- count_matrix[-idx,]

## Sort annotations to match order in counts
annotation_df %<>%
  arrange(
    match(specimenID, colnames(count_matrix))
  )

# Make sure that counts and annotations have the same order
is_identical <- rownames(annotation_df) == colnames(count_matrix)
if (!all(is_identical)) {
  stop("There are unannotated samples in the counts", call. = FALSE)
}

# Change numeric codes to factors, add explanatory
# level names to some, rename columns, and remove unneeded columns

## Use only rownames (specimenID) as keys
annotation_df %<>%
  mutate(
    specimenID = NULL,
    projid = NULL,
    individualID = NULL,
    ID = NULL,
    Sampleid = NULL
  )

## Shorten picard data column names (PF_ALIGNED_BASES is a special case),
## perl is required for regex lookahead
colnames(annotation_df) %<>%
  gsub('AlignmentSummaryMetrics__(?!PF_ALIGNED_BASES)', '', ., perl = TRUE) %>%
  gsub('AlignmentSummaryMetrics__', 'Alignment_', .) %>%
  gsub('RnaSeqMetrics__(?!PF_ALIGNED_BASES)', '', ., perl = TRUE) %>%
  gsub('RnaSeqMetrics__', 'RnaSeq_', .)

## Combine information from libraryBatch and Batch
annotation_df$libraryBatch %<>% as.factor
annotation_df$Batch %<>%
  as.factor %>%
  coalesce(annotation_df$libraryBatch) %>%
  droplevels

annotation_df %<>%
  mutate(libraryBatch = NULL)

annotation_df$Study %<>% as.factor
annotation_df$msex %<>% as.factor
levels(annotation_df$msex) <- c("Female", "Male")

annotation_df$race %<>% as.factor
annotation_df$spanish %<>% as.factor
annotation_df$apoe_genotype %<>% as.factor
annotation_df$braaksc %<>% as.factor

annotation_df$ceradsc %<>% as.factor
levels(annotation_df$ceradsc) <- c("Definite", "Probable", "Possible", "No_AD")

annotation_df$cogdx %<>% as.factor
levels(annotation_df$cogdx) <- c("NCI", "MCI", "MCI+", "AD", "AD+", "Other")

annotation_df$dcfdx_lv %<>% as.factor
levels(annotation_df$dcfdx_lv) <- c("NCI", "MCI", "MCI+", "AD", "AD+", "Other")

## Add a binary ceradsc
annotation_df$ceradsc_binary <- annotation_df$ceradsc
levels(annotation_df$ceradsc_binary) <- c("AD", "AD", "No_AD", "No_AD")

## Change age from string to numeric
annotation_df$age_at_visit_max %<>%
  sub("\\+", "", .) %>% # 90+ is treated as 90
  as.numeric
annotation_df$age_death %<>%
  sub("\\+", "", .) %>% # 90+ is treated as 90
  as.numeric

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
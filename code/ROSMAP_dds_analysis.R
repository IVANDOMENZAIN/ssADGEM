#!/usr/bin/env Rscript

# TODO: Experiment with cqn normalisation
# TODO: Remove outliers not in batch 7,8
# TODO: Rigorously do covariate significance testing

# Parse arguments from terminal
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Format: ./<script> <dds-file> <target>", call. = FALSE)
  }
  dds_path <- args[1]
  target_folder <- args[2]
} else {
  # For using interactively, change paths to appropriate
  project_path <- "/castor/project/home/gasto/ssADGEM/"
  dds_path <- paste(project_path, "data/ROSMAP_dds.rds", sep = "")
  target_folder <- paste(project_path, "data/", sep = "")
}


# Needed libraries
library("BiocParallel")
register(MulticoreParam())
library("tidyverse")
library("magrittr") # Piping
library("cqn") # Normalise for future use
if (interactive()){
  library("DESeq2")
  library("ggplot2")
  library("ggforce")
  library("pheatmap")
  library("limma") # removeBatchEffect
  library("broom") # PCA
}

# Load the dataset
load(dds_path)

# Read gene parameters
working_directory = getwd()
"external/synapse_dir/" %>%
  paste(project_path, ., sep = "") %>%
  setwd

gene_data <- "geneParameters.tsv" %>%
  read_delim(delim = "\t") %>%
  as.data.frame

gene_data %<>%
  mutate(hgnc_symbol = NULL) %>% # No way to know which one is meant in counts
  filter(!is.na(gene.length)) %>%
  filter(!is.na(percentage_gc_content)) %>%
  distinct

working_directory %>% setwd

# Exclude samples with other cognitive conditions
idx_other_CI <- which(
  dds$cogdx == "MCI+" |
    dds$cogdx == "AD+" |
    dds$cogdx == "Other")
dds <- dds[, -idx_other_CI]
dds$cogdx %<>% droplevels

# Analysis (only if interactive)
if (interactive()){
  
  # Pre-filter out genes with CPM < 1 in more than 
  # half of the samples for a given diagnosis
  genes_to_keep <- c()
  for (dx in c("NCI", "MCI", "AD")) {
    
    # Get the respective counts
    cpm <- dds %>% .[, .$cogdx == dx] %>%
      counts
    
    # Convert to CPM
    cpm[,] %<>% scale(center = FALSE, scale = colSums(.))
    cpm <- cpm * 1e6
    
    # Filter
    genes_non_low_percentage <- rowMeans(cpm > 1)
    genes_to_keep <- rownames(cpm)[genes_non_low_percentage > .5] %>%
      union(genes_to_keep)
    
  }
  genes_to_keep %<>% # Exclude genes without gene data
    intersect(gene_data$Gene.ID)
  gene_data.filtered <- gene_data %>% filter(Gene.ID %in% genes_to_keep)
  dds.filtered <- dds[genes_to_keep, ]
  
  # Cross-tabulation of cogdx and CERAD
  
  comp_df <- dds %>%
    colData %>%
    as.data.frame
  comp_df %<>% .[, c("cogdx", "ceradsc")]
  comp_df %<>%
    table %>%
    prop.table(margin = "cogdx") %>%
    as.data.frame(responseName = "percentcogdx")
  
  comp_df %>% ggplot(aes(x = ceradsc, y = cogdx)) +
    geom_tile(aes(fill = percentcogdx, height = .9)) +
    geom_text(aes(label = percentcogdx %>% round(2)), colour = "white") +
    scale_fill_continuous("cogdx%")
  
  # Show that batch is a significant co-variate from PCA plot
  
  cqnAndDESeq <- function(deseqds, gene_param, dsgn) {
    # Set design
    design(deseqds) <- dsgn
    
    # Estimate size factors with Conditional Quantile
    stopifnot(rownames(deseqds) == gene_param$Gene.ID)
    cqn_res <- deseqds %>% assay %>%
      cqn(x = gene_param$percentage_gc_content,
          lengths = gene_param$gene.length)
    normalizationFactors(deseqds) <- cqn_res$glm.offset %>% exp
    
    # Perform rest of DeSeq pipeline
    deseqds %<>% estimateDispersions(parallel = TRUE)
    deseqds %<>% nbinomWaldTest(parallel = TRUE)
    
    return(deseqds)
  }
  
  ## Set design as AD or not and perform estimations
  dds.filtered %<>% cqnAndDESeq(gene_data.filtered, ~ cogdx)
  design(dds.filtered) <- ~ cogdx
  dds.filtered %<>% DESeq(parallel = TRUE)
  
  ## Transform (not with rlog due to speed)
  res <- dds.filtered %>% results
  vsd <- dds.filtered[res$padj < .1, ] %>%
    vst(blind = FALSE)
  
  ## Show clustering by batch of significant genes
  pcaFromVSD <- function(dst, var = c("Batch", "msex")) {
    
    count_df <- dst %>%
      assay %>%
      t %>%
      as.data.frame
    
    metadata <- dst %>% colData %>% as.data.frame
    
    count_df %<>% merge(metadata,
                        by = "row.names")
    
    pca_fit <- count_df %>%
      dplyr::select(starts_with("ENSG")) %>%
      prcomp
    
    pcPercent <- pca_fit %>%
      tidy(matrix = "eigenvalues") %>%
      .$percent
    
    plt <- pca_fit %>%
      augment(count_df) %>%
      ggplot(aes(.fittedPC1, .fittedPC2, color = get(var[1]), shape = get(var[2]))) +
      geom_point(size = 3) +
      labs(color = var[1], shape = var[2]) +
      xlab(sprintf("PC1: %2.0f%% variance",100*pcPercent[1])) +
      ylab(sprintf("PC2: %2.0f%% variance",100*pcPercent[2]))
    return(plt)
  }
  pcaFromVSD(vsd)
  
  # Because batch 7 and 8 have an outsized effect on the data, we remove them
  dds.sans78 <- dds.filtered[, which(dds.filtered$Batch != 7 & dds.filtered$Batch != 8)]
  dds.sans78$Batch %<>% droplevels
  
  #Estimate new parameters
  design(dds.sans78) <- ~ cogdx
  dds.sans78 %<>% DESeq(parallel = TRUE)
  
  # Compare new PCA and heatmap, we see some batch effect
  res <- dds.sans78 %>% results
  vsd.sans78 <- dds.sans78[res$padj < .1, ] %>%
    vst(blind = FALSE)
  pcaFromVSD(vsd.sans78)
  
  ## Show batch and RIN clustering on heatmap
  heatmapFromVSD <- function(dsds, dst, var = c("cogdx", "ceradsc_binary")) {
    
    res <- dsds %>%
      results(contrast = c("cogdx", "AD", "NCI"))
    
    idx_signGenes <- which(res$padj < .1)
    dst <- dst[idx_signGenes, ]
    
    normCounts <- dst %>% assay
    nGenes <- 1000
    if (length(idx_signGenes) < nGenes)
      nGenes <- length(idx_signGenes)
    title <- sprintf("Heatmap of top %s variable DE genes (padj < .1)", nGenes)
    
    topVarGenes <- order(-rowVars(normCounts)) %>% head(nGenes)
    
    mat <- assay(dst)[topVarGenes, ]
    mat <- mat - rowMeans(mat)
    
    df <- as.data.frame(colData(dst)[, var])
    pheatmap(mat,
             main = title,
             annotation_col = df,
             show_colnames = FALSE,
             show_rownames = FALSE)
  }
  heatmapFromVSD(dds.sans78[res$padj < .1, ], vsd.sans78, c("RIN", "Batch"))
  pcaFromVSD(vsd.sans78, c("RIN", "msex"))
  
  ## The RIN <-> batch correlation is reflected in the methods used
  ### According to the article responsible for RNA-seq (doi:10.1038/s41593-018-0154-9):
  ### "The libraries were constructed and pooled according to the RIN
  ### scores such that similar RIN scores would be pooled together.
  ### Varying RIN scores results in a larger spread of insert sizes
  ### during library construction and leads to uneven coverage distribution
  ### throughout the pool."
  
  # Thus we include RIN and batch as covariates
  dds.covar <- dds.sans78
  design(dds.covar) <- ~ Batch + RIN + cogdx
  dds.covar %<>% DESeq(parallel = TRUE)
  
  # Show the effect of correcting for batch, msex and RIN
  vsd.covar <- dds.covar[res$padj < .1, ] %>%
    vst(blind = FALSE)
  pcaFromVSD(vsd.covar, c("RIN", "msex"))
  heatmapFromVSD(dds.covar[res$padj < .1, ], vsd.covar, c("RIN", "Batch"))
  
  batchAndRINCorrect <- function(dst, dsgn) {
    mm <- model.matrix(dsgn, data = colData(dst))
    assay(dst) %<>% removeBatchEffect(
      batch = dst$libraryBatch,
      covariates = dst$RIN,
      design = mm)
    return(dst)
  }
  design(dds.sans78) <- ~ Batch + RIN + cogdx
  dds.sans78 %<>% DESeq(parallel = TRUE)
  vsd.corrected <- dds.sans78 %>%
    vst(blind = FALSE) %>%
    batchAndRINCorrect(~ cogdx)
  pcaFromVSD(vsd.corrected, var = c("Batch", "msex"))
  heatmapFromVSD(dds.sans78, vsd.corrected)
  
} else { # Not interactive
  
  # Normalisation
  
  ## Subset to only supported
  dds.supported <- dds[gene_data$Gene.ID,
                       which(dds$Batch != 7 & dds$Batch != 8)]
  dds.supported$Batch %<>% droplevels
  
  (rownames(dds.supported) == gene_data$Gene.ID) %>% stopifnot
  
  ## Run DeSeq2
  design(dds.supported) <- ~ cogdx
  dds.supported %<>% DESeq(parallel = TRUE)
  
  ## Normalise by previous size factors
  vsd.supported <- dds.supported %>%
    vst(blind = FALSE)

  ## Correct for batch and RIN
  mat <- vsd.supported %>% assay
  mm <- model.matrix(~ cogdx,
                     vsd.supported %>% colData)
  mat <- removeBatchEffect(mat,
                           batch = vsd.supported$Batch,
                           covariates = vsd.supported$RIN,
                           design = mm)
  
  ## Export
  target_path_counts <- target_folder %>%
    paste("ROSMAP_normalized_log2counts.txt", sep = "")
  mat %>% write.table(file = target_path_counts, col.names = NA)

  target_path_meta <- target_folder %>% paste("ROSMAP_metadata.txt", sep = "")
  dds.supported %>% colData %>%
    write.table(file = target_path_meta, col.names = NA)
  
}
# Caret for copy+paste: ^
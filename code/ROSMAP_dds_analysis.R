#!/usr/bin/env Rscript

# TODO: Remove outliers not in batch 7
# TODO: Rigorously do covariate significance testing


 ############
 # PREAMBLE #
 ############


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
  plt_folder <- paste(project_path, "nobackup/plots/", sep = "")
  target_folder <- paste(project_path, "data/", sep = "")
}


# Needed libraries
library("BiocParallel")
register(MulticoreParam())
library("tidyverse")
library("magrittr") # Piping
library("cqn")
if (interactive()){
  library("DESeq2")
  library("ggplot2")
  library("ggforce")
  library("pheatmap")
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


 #################
 # DATA ANALYSIS #
 #################


# Exclude samples with other cognitive conditions
idx_other_CI <- which(
  dds$cogdx == "MCI+" |
  dds$cogdx == "AD+" |
  dds$cogdx == "Other")
dds <- dds[, -idx_other_CI]
dds$cogdx %<>% droplevels

# Pre-filter out genes with counts below 10 in at least minCohortSize samples
minCohortSize <- dds$cogdx %>%
  table %>%
  min
genes_to_keep <- rowSums(assay(dds) >= 100) >= minCohortSize
dds.filtered <- dds[genes_to_keep, ]

# Also filter out genes without gene data
gene_data.filtered <- gene_data %>% filter(Gene.ID %in% rownames(dds.filtered))
dds.filtered <- dds[gene_data.filtered$Gene.ID, ]

# Analysis (only if interactive)
if (interactive()){

  
  ###################
  # DEA STARTS HERE #
  ###################
  
  # Helper function
  plot.new()
  savePlot <- function(new_plt, save_path = "plt%03d.jpeg") {
    jpeg(filename = save_path,
         pointsize = 12,
         height = 600,
         width = 900,
         quality = 95)
    plot(new_plt)
    invisible(dev.off())
  }
  
  # Check the distributions (logCPM)
  checkDistribution <- function(data, plt_title = "Gene distributions per sample"){
    
    if (is(data, "DESeqDataSet")){
      cpm <- data %>%
        fpm(robust = FALSE)
      logCPM <- log2(0.01+cpm) # avoid -Inf
      x_title <- bquote(log[2](0.01+CPM))
    } else {
      logCPM <- data
      x_title <- bquote(log[2](CPM))
    }
    plt <- logCPM %>%
      as.data.frame %>%
      stack %>%
      ggplot(aes(x = values, color = ind)) +
      geom_density(aes(group = ind), linetype = "dashed") +
      ggtitle(plt_title) +
      xlab(x_title) +
      theme(legend.position = "none")
    return(plt)
  }
  plt_path <- paste(plt_folder, "gene_dist_high_filter.jpeg", sep = "")
  dds.filtered %>%
    checkDistribution(plt_title = "Gene distribution of all samples (high pre-filtration)") %>%
    savePlot(save_path = plt_path)
  
  # Distributions do not have good overlap, use cqn to normalise
  cqn_res <- dds.filtered %>% assay %>%
    cqn(x = gene_data.filtered$percentage_gc_content,
        lengths = gene_data.filtered$gene.length)
  
  plt_path <- paste(plt_folder, "gene_dist_cqn_high_filter.jpeg", sep = "")
  (cqn_res$y + cqn_res$offset) %>%
    checkDistribution(plt_title = "Gene distribution after cqn (high pre-filtration)") %>%
    savePlot(save_path = plt_path)
  
  cqnNormalizationFactors <- function(deseqds, gene_param) {
    # Estimate size factors with Conditional Quantile
    stopifnot(rownames(deseqds) == gene_param$Gene.ID)
    cqn_res <- deseqds %>% assay %>%
      cqn(x = gene_param$percentage_gc_content,
          lengths = gene_param$gene.length)
    normalizationFactors(deseqds) <- cqn_res$glm.offset %>% exp
    return(deseqds)
  }
  getVSDFromDeSeqDS <- function(deseqds, dsgn){
    design(deseqds) <- dsgn
    dsds %<>% cqnNormalizationFactors(gene_data.filtered)
    dsds %<>% DESeq(parallel = TRUE)
    res <- dsds %>% results
    vsd <- dsds[res$padj < .1, ] %>%
      vst(blind = FALSE)
    return(vsd)
  }
  
  ## Set design as AD or not and perform estimations
  dds.filtered %<>% .[gene_data.filtered$Gene.ID, ]
  vsd.filtered <- dds.filtered %>%
    getVSDFromDeSeqDS(~ cogdx)
  
  
  ######################################
  # SHOWING THAT BATCH 7 IS AN OUTLIER #
  ######################################
  
  
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
  heatmapFromVSD <- function(dst, var = c("cogdx", "Batch")) {
    
    normCounts <- dst %>% assay
    nGenes <- min(normCounts %>% dim %>% .[1],
                  1000)
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
  vsd.filtered %>% pcaFromVSD
  vsd.filtered %>% heatmapFromVSD

  ## Because batch 7 has an outsized effect on the data, we remove it
  dds.sans7 <- dds.filtered[, which(dds.filtered$Batch != 7)]
  dds.sans7$Batch %<>% droplevels

  vsd.sans7 <- dds.sans7 %>%
    getVSDFromDeSeqDS(~ cogdx)

  
  ###################################
  # SHOWING THAT RIN IS A COVARIATE #
  ###################################
  
  
  vsd.sans7 %>% pcaFromVSD(var = c("RIN", "msex"))
  vsd.sans7 %>% heatmapFromVSD(var = c("Batch", "RIN"))
  
  ## The RIN <-> batch correlation is reflected in the methods used
  ### According to the article responsible for RNA-seq (doi:10.1038/s41593-018-0154-9):
  ### "The libraries were constructed and pooled according to the RIN
  ### scores such that similar RIN scores would be pooled together.
  ### Varying RIN scores results in a larger spread of insert sizes
  ### during library construction and leads to uneven coverage distribution
  ### throughout the pool."
  
  # Thus we include RIN and batch as covariates. This will not affect
  # normalisation, but will affect which genes are considered significant and
  # therefore which genes are used to cluster the genes
  dds.covar <- dds.sans7
  vsd.covar <- dds.covar %>%
    getVSDFromDeSeqDS(~ Batch + RIN + cogdx)
  
  vsd.covar %>% pcaFromVSD(var = c("RIN", "msex"))
  vsd.covar %>% heatmapFromVSD(var = c("RIN", "Batch"))
  
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
    paste("ROSMAP_normalized_log2counts.txt.gz", sep = "")
  mat %>% write.table(file = gzfile(target_path_counts), col.names = NA)
  
  target_path_counts <- target_folder %>% paste("ROSMAP_metadata.txt.gz", sep = "")
  dds.supported %>% colData %>%
    write.table(file = gzfile(target_path_counts), col.names = NA)
  
}
# Caret for copy+paste on thinlinc: ^
# Tilde for copy+paste on thinlinc: ~

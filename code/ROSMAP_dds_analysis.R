#!/usr/bin/env Rscript

# TODO: Remove outliers by SD
# TODO: Do covariate significance testing in a rigorous manner

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
  library("latentcor") # Correlation Plot
}


 ##############
 # LOAD FILES #
 ##############


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
genes_to_keep <- rowSums(assay(dds) >= 10) >= minCohortSize
dds.filtered <- dds[genes_to_keep, ]

# Also filter out genes without gene data
gene_data.filtered <- gene_data %>% filter(Gene.ID %in% rownames(dds.filtered))
dds.filtered <- dds[gene_data.filtered$Gene.ID, ]

# Add rowRanges
tmp <- dds.filtered %>% rownames
gr <- GRanges(seqnames = gene_data.filtered$Gene.ID %>% as.character,
        ranges = IRanges(start = 1, width = gene_data.filtered$gene.length))
rowRanges(dds.filtered) <- gr %>% split
rownames(dds.filtered) <- tmp


 ###################
 # DEA STARTS HERE #
 ###################

# Helper functions
plot.new()
savePlot <- function(new_plt, save_path = "plt%03d.jpeg") {
  if(new_plt %>% is.list) {
    jpeg(filename = save_path,
         pointsize = 12,
         height = 600,
         width = 900,
         quality = 95)
    grid.draw(new_plt$gtable)
    invisible(dev.off())
  } else {
    jpeg(filename = save_path,
         pointsize = 12,
         height = 600,
         width = 900,
         quality = 95)
    plot(new_plt)
    invisible(dev.off())
  }
}
checkpoint <- function (obj_name) {
  project_path <- "/castor/project/home/gasto/ssADGEM/"
  checkpoint_path <- paste(project_path, "nobackup/chpt_",
                           obj_name, ".Rdata", sep = "")
  if (exists(obj_name)) {
    save.image(file = checkpoint_path)
  } else {
    load(file = checkpoint_path, envir = globalenv())
  }
}

if (interactive()){

  ########################
  # METADATA CORRELATION #
  ########################
  
  # Calculate correlations
  metadata_estimates <- dds.filtered %>% colData %>%
    latentcor(types = get_types(.))
  
  checkpoint("metadata_estimates")
  
  # Plot pointwise correlation
  plt_path <- paste(plt_folder, "metadata_assoc.jpeg", sep = "")
  metadata_estimates$Rpointwise %>% abs %>%
    pheatmap(main = "Metadata Association Heatmap",
             color = viridis::viridis(8),
             filename = plt_path)

  ########################
  # PER BATCH RLE AND QN #
  ########################
  
  batches <- dds.filtered$Batch %>% levels
  
  # Plotting functions
  rle_plot <- function(expression_data,
                       TRANSFORM = function(x) log2(1+x)
                       ) {
    
    log_expression <- TRANSFORM(expression_data)
    
    gene_log_median <- log_expression %>%
      apply(MARGIN = 1, median)
    
    relative_log_expression <- log_expression %>%
      sweep(MARGIN = 1, STATS = gene_log_median)
    
    relative_log_expression %>%
      as.data.frame %>%
      pivot_longer(cols = everything(),
                   names_to = "sample",
                   values_to = "rle") %>%
      ggplot(aes(x = sample, y = rle)) +
      geom_boxplot(fill = "slateblue", alpha = 0.2, outlier.shape = NA) +
      coord_cartesian(ylim = c(-5,5)) +
      theme(axis.text.x = element_blank())%>%
      return()
  }
  
  # RLE plot per batch before quantile normalisation
  for (batch in batches) {
    
    plt_path <- paste(plt_folder, "raw_rle_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Relative Log Expression in batch", batch,
                   "before quantile normalisation")
    
    idx <- dds.filtered$Batch == batch
    dds.filtered[, idx] %>% assay %>%
      rle_plot(TRANSFORM = function(x) log2(1+x)) %>%
      savePlot(save_path = plt_path)
    
  }
  
  # Quantile Normalisation
  perBatchQuantile <- function(deseqds) {
    deseqds.new <- deseqds
    for (batch in batches) {
      idx <- deseqds.new$Batch == batch
      assay(deseqds.new[, idx], withDimnames = FALSE) %<>%
        normalize.quantiles(keep.names = TRUE)
    }
    return(deseqds.new)
  }
  dds.qn <- dds.filtered %>% perBatchQuantile
  
  # RLE plot per batch after quantile normalisation
  for (batch in batches) {
    
    plt_path <- paste(plt_folder, "post_qn_rle_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Relative Log Expression in batch", batch,
                   "after quantile normalisation")
    
    idx <- dds.qn$Batch == batch
    dds.qn[, idx] %>% assay %>%
      rle_plot(TRANSFORM = function(x) log2(1+x)) %>%
      savePlot(save_path = plt_path)
    
  }
  
  ######################
  # APPLY FPKM AND TPM #
  ######################
  
  checkDistribution <- function(data,
                                plt_title = "Gene distributions per sample",
                                xlabel = bquote(log[2](1+CPM))
  ) {
    
    if (is(data, "DESeqDataSet")){
      cpm <- data %>%
        fpm(robust = FALSE)
      logCPM <- log2(1+cpm)
    } else {
      logCPM <- data
    }
    plt <- logCPM %>%
      as.data.frame %>%
      stack %>%
      ggplot(aes(x = values, color = ind)) +
      geom_density(aes(group = ind), linetype = "dashed") +
      ggtitle(plt_title) +
      xlab(xlabel) +
      theme(legend.position = "none")
    return(plt)
  }
  
  # Pre- versus post-qn CPM distributions per batch
  for (batch in batches) {
    
    idx <- dds.qn$Batch == batch
    
    # Pre
    plt_path <- paste(plt_folder, "gene_dist_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution in batch", batch,
                   "before quantile normalisation")
    
    dds.filtered[, idx] %>%
      checkDistribution(plt_title = title) %>%
      savePlot(save_path = plt_path)
    
    # Post
    plt_path <- paste(plt_folder, "gene_dist_post_qn_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution in batch", batch,
                   "after quantile normalisation")
    
    dds.qn[, idx] %>%
      checkDistribution(plt_title = title) %>%
      savePlot(save_path = plt_path)
    
  }
  
  # FPKM and TPM + Quantile Normalisation
  dds.fpkm <- dds.filtered
  assay(dds.fpkm, withDimnames = FALSE) <- dds.filtered %>%
    fpkm(robust = FALSE)
  
  dds.tpm <- dds.fpkm
  assay(dds.tpm) <- assay(dds.fpkm) %>%
    sweep(2, colSums(.), '/')
  assay(dds.tpm) <- 1e6 * assay(dds.tpm)
  
  dds.fpkm %<>% perBatchQuantile
  dds.tpm %<>% perBatchQuantile
  
  for (batch in batches) {
    
    idx <- dds.filtered$Batch == batch
    
    # FPKM
    plt_path <- paste(plt_folder, "gene_dist_fpkm_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution in batch", batch,
                   "after transforming to FPKM")
    log2(1+assay(dds.fpkm)[, idx]) %>%
      checkDistribution(plt_title = title,
                        xlabel = bquote(log[2](1+FPKM))) %>%
      savePlot(save_path = plt_path)
    
    # TPM
    plt_path <- paste(plt_folder, "gene_dist_tpm_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution in batch", batch,
                   "after transforming to TPM")
    
    log2(1+assay(dds.tpm)[, idx]) %>%
      checkDistribution(plt_title = title,
                        xlabel = bquote(log[2](1+TPM))) %>%
      savePlot(save_path = plt_path)
    
  }
  
  #############
  # RUN DESEQ #
  #############

  runDESeqWithDesign <- function(deseqds, dsgn){
    design(deseqds) <- dsgn
    deseqds %>%
      DESeq(parallel = TRUE) %>%
      return()
      
  }
  getVSDByPadj <- function(deseqds, min_padj = .1) {
    res <- deseqds %>% results
    deseqds[res$padj < min_padj, ] %>%
      vst(blind = FALSE) %>%
      return()
  }
  
  # Set design as AD or not and perform estimations
  dds.filtered %<>% runDESeqWithDesign(~ cogdx)
  
  # Get variance stabilized values
  vsd.filtered <- dds.filtered %>%
    getVSDByPadj(min_padj = 0.1)
  
  #################################
  # SHOW BATCH EFFECTS OF 7 AND 8 #
  #################################
  
  checkpoint("vsd.filtered")
  
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
    plt <- pheatmap(mat,
             main = title,
             annotation_col = df,
             show_colnames = FALSE,
             show_rownames = FALSE)
    return(plt)
  }
  
  plt_path <- paste(plt_folder, "significant_genes_pca_plot.jpeg", sep = "")
  vsd.filtered %>%
    pcaFromVSD %>% savePlot(save_path = plt_path)
  
  plt_path <- paste(plt_folder, "significant_genes_heatmap.jpeg", sep = "")
  a <- vsd.filtered %>%
    heatmapFromVSD %>% savePlot(save_path = plt_path)
  
  # Set new design including batch as covariate
  dds.batch_design <- dds.filtered %>% runDESeqWithDesign(~ Batch + cogdx)
  vsd.batch_design <- dds.batch_design %>%
    getVSDByPadj(min_padj = 0.1)
  
  ###################################
  # SHOWING THAT RIN IS A COVARIATE #
  ###################################
  
  checkpoint("vsd.batch_design")
  
  plt_path <- paste(plt_folder, "rin_pca_plot.jpeg", sep = "")
  vsd.batch_design %>%
    pcaFromVSD(var = c("RIN", "msex")) %>% savePlot(save_path = plt_path)
  
  plt_path <- paste(plt_folder, "rin_heatmap.jpeg", sep = "")
  vsd.batch_design %>%
    heatmapFromVSD(var = c("Batch", "RIN")) %>% savePlot(save_path = plt_path)
  
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
  # Set new design including batch as covariate
  dds.rin_design <- dds.filtered %>% runDESeqWithDesign(~ RIN + Batch + cogdx)
  vsd.rin_design <- dds.rin_design %>%
    getVSDByPadj(min_padj = 0.1)
  
  checkpoint("vsd.rin_design")
  
  plt_path <- paste(plt_folder, "covar_design_pca_plot.jpeg", sep = "")
  vsd.rin_design %>%
    pcaFromVSD(var = c("RIN", "msex")) %>% savePlot(save_path = plt_path)
  
  plt_path <- paste(plt_folder, "covar_design_heatmap.jpeg", sep = "")
  vsd.rin_design %>%
    heatmapFromVSD(var = c("Batch", "RIN")) %>% savePlot(save_path = plt_path)
  
} else {
  #############################################
  # NOT INTERACTIVE, EXPORT NORMALISED COUNTS #
  #############################################
  
  runDESeqWithDesign <- function(deseqds, dsgn){
    design(deseqds) <- dsgn
    deseqds %>%
      DESeq(parallel = TRUE) %>%
      return()
    
  }
  getVSDByPadj <- function(deseqds, min_padj = .1) {
    res <- deseqds %>% results
    deseqds[res$padj < min_padj, ] %>%
      vst(blind = FALSE) %>%
      return()
  }
  
  # Set design as AD or not and perform estimations
  dds.filtered %<>% runDESeqWithDesign(~ cogdx)
  
  # Get variance stabilized values
  vsd.filtered.all_genes <- dds.filtered %>%
    getVSDByPadj(min_padj = 1)
  
  # Export
  target_path_counts <- target_folder %>%
    paste("ROSMAP_normalized_log2counts.txt.gz", sep = "")
  vsd.filtered.all_genes %>%
    assay %>%
    write.table(file = gzfile(target_path_counts),
                col.names = NA)
  
  target_path_meta <- target_folder %>%
    paste("ROSMAP_metadata.csv.gz", sep = "")
  dds.filtered %>% colData %>%
    write.table(file = gzfile(target_path_meta),
                col.names = NA)
  
}
# Caret for copy+paste on thinlinc: ^
# Tilde for copy+paste on thinlinc: ~

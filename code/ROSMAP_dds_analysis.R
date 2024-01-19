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
  library("grid")
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

# Plot save function
plot.new()
savePlot <- function(new_plt, save_path = "plt%03d.jpeg") {
  
  # Save as small
  jpeg(filename = save_path,
       pointsize = 12,
       height = 300,
       width = 450,
       quality = 95)
  plot(new_plt)
  invisible(dev.off())
  
  save_path %<>% 
    strsplit(split = ".jpeg") %>%
    paste("_large.jpeg", sep = "")
  
  # Save as large
  jpeg(filename = save_path,
       pointsize = 12,
       height = 600,
       width = 900,
       quality = 95)
  plot(new_plt)
  invisible(dev.off())
}

if (interactive()){

  ########################
  # METADATA CORRELATION #
  ########################
  
  # Calculate correlations
  metadata_estimates <- dds.filtered %>% colData %>%
    latentcor(types = get_types(.))
  
  # Plot pointwise correlation
  plt_path <- paste(plt_folder, "metadata_assoc.jpeg", sep = "")
  metadata_estimates$Rpointwise %>% abs %>%
    pheatmap(main = "Metadata Association Heatmap",
             color = viridis::viridis(8),
             filename = plt_path,
             width = 9,
             height = 6)

  ########################
  # PER BATCH RLE AND QN #
  ########################
  
  batches <- dds.filtered$Batch %>% levels
  batch_key <- dds.filtered %>%
    colData %>%
    as.data.frame %>%
    select(Batch)
  cogdx_key <- dds.filtered %>%
    colData %>%
    as.data.frame %>%
    select(cogdx)
  
  # Plotting functions
  rle_plot <- function(expression_data,
                       TRANSFORM = function(x) log2(1+x),
                       plt_title = ""
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
      mutate(cogdx = cogdx_key[sample, ]) %>%
      ggplot(aes(x = sample, y = rle, color = cogdx)) +
      geom_boxplot(alpha = 0.2, outlier.shape = NA) +
      labs(title = plt_title) +
      coord_cartesian(ylim = c(-5,5)) +
      theme(axis.text.x = element_blank(), legend.position = "none") +
      facet_wrap(vars(cogdx), scales = "free_x") %>%
      return()
  }
  rle_plot_by_batch <- function(expression_data,
                       TRANSFORM = function(x) log2(1+x),
                       plt_title = ""
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
      mutate(batch = batch_key[sample, ]) %>%
      ggplot(aes(x = sample, y = rle, color = batch)) +
      geom_boxplot(alpha = 0.2, outlier.shape = NA) +
      labs(title = plt_title) +
      coord_cartesian(ylim = c(-5,5)) +
      theme(axis.text.x = element_blank(), legend.position = "none") +
      facet_wrap(vars(batch), scales = "free_x") %>%
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
      rle_plot(TRANSFORM = function(x) log2(1+x),
               plt_title = title) %>%
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
  
  plt_path <- paste(plt_folder, "post_qn_rle_batch.jpeg", sep = "")
  title <- paste("Relative Log Expression in after quantile normalisation")
  dds.qn %>% assay %>%
    rle_plot_by_batch(TRANSFORM = function(x) log2(1+x), plt_title = title) %>%
    savePlot(save_path = plt_path)
  
  ##################
  # APPLY TPM + QN #
  ##################
  
  checkDistribution <- function(data,
                                plt_title = "Gene distributions per sample",
                                xlabel = bquote(log[2](1+TPM))
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

  # TPM
  dds.fpkm <- dds.filtered
  assay(dds.fpkm, withDimnames = FALSE) <- dds.filtered %>%
    fpkm(robust = FALSE)
  
  dds.tpm <- dds.fpkm
  assay(dds.tpm) <- assay(dds.fpkm) %>%
    sweep(2, colSums(.), '/')
  assay(dds.tpm) <- 1e6 * assay(dds.tpm)
  
  for (batch in batches) {
    
    idx <- dds.filtered$Batch == batch
    
    # Pre
    plt_path <- paste(plt_folder, "gene_dist_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution in batch", batch,
                   "before TPM")
    dds.filtered[, idx] %>%
      checkDistribution(plt_title = title) %>%
      savePlot(save_path = plt_path)
    
    # Post
    plt_path <- paste(plt_folder, "gene_dist_TPM_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution in batch", batch,
                   "after TPM")
    dds.tpm[, idx] %>%
      checkDistribution(plt_title = title) %>%
      savePlot(save_path = plt_path)
    
  }
  
  # Quantile Normalisation
  dds.tpm.qn <- dds.tpm %>% perBatchQuantile
  
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
    res$padj[res$padj %>% is.na] <- 1
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
      geom_point(size = 2, alpha = .75) +
      labs(color = var[1], shape = var[2]) +
      xlab(sprintf("PC1: %2.0f%% variance",100*pcPercent[1])) +
      ylab(sprintf("PC2: %2.0f%% variance",100*pcPercent[2]))
    return(plt)
  }
  
  plt_path <- paste(plt_folder, "significant_genes_batch_pca_plot.jpeg", sep = "")
  vsd.filtered %>%
    pcaFromVSD %>% savePlot(save_path = plt_path)
  
  # Set new design including batch as covariate
  dds.batch_design <- dds.filtered %>% runDESeqWithDesign(~ Batch + cogdx)
  vsd.batch_design <- dds.batch_design %>%
    getVSDByPadj(min_padj = 0.1)
  
  ###################################
  # SHOWING THAT RIN IS A COVARIATE #
  ###################################
  
  plt_path <- paste(plt_folder, "significant_genes_rin_pca_plot.jpeg", sep = "")
  vsd.batch_design %>%
    pcaFromVSD(var = c("RIN", "ceradsc_binary")) %>%
    savePlot(save_path = plt_path)
  
  ## The RIN <-> batch correlation is reflected in the methods used
  ### According to the article responsible for RNA-seq (doi:10.1038/s41593-018-0154-9):
  ### "The libraries were constructed and pooled according to the RIN
  ### scores such that similar RIN scores would be pooled together.
  ### Varying RIN scores results in a larger spread of insert sizes
  ### during library construction and leads to uneven coverage distribution
  ### throughout the pool."
  
  # Plot RIN-Batch relation
  plt_path <- paste(plt_folder, "rin2batch_boxplot.jpeg", sep = "")
  plt <- dds.filtered %>%
    colData %>%
    as.data.frame %>%
    ggplot(aes(x=Batch, y=RIN, fill=Batch)) +
    geom_boxplot() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("RIN distribution by batch")
  plt %>%
    savePlot(save_path = plt_path)

  # Thus we include RIN and batch as covariates. This will not affect
  # normalisation, but will affect which genes are considered significant and
  # therefore which genes are used to cluster the genes
  # Set new design including batch as covariate
  dds.rin_design <- dds.filtered %>% runDESeqWithDesign(~ RIN + Batch + cogdx)
  vsd.rin_design <- dds.rin_design %>%
    getVSDByPadj(min_padj = 0.1)
  
  
  ##################################
  # SEEING IF AD DIAGNOSES CLUSTER #
  ##################################
  
  # Check cogdx and ceradsc
  plt_path <- paste(plt_folder,
                    "significant_genes_cogdx_pca_plot.jpeg", sep = "")
  vsd.rin_design %>%
    pcaFromVSD(var = c("cogdx", "ceradsc")) %>% savePlot(save_path = plt_path)

  # Check ceradsc_binary
  plt_path <- paste(plt_folder,
                    "significant_genes_ceradsc_binary_pca_plot.jpeg", sep = "")
  vsd.rin_design %>%
    pcaFromVSD(var = c("ceradsc_binary", "ceradsc_binary")) %>%
    savePlot(save_path = plt_path)

} else {
  #############################################
  # NOT INTERACTIVE, EXPORT NORMALISED COUNTS #
  #############################################
  
  batches <- dds.filtered$Batch %>% levels
  
  # TPM + Quantile Normalisation
  perBatchQuantile <- function(deseqds) {
    deseqds.new <- deseqds
    for (batch in batches) {
      idx <- deseqds.new$Batch == batch
      assay(deseqds.new[, idx], withDimnames = FALSE) %<>%
        normalize.quantiles(keep.names = TRUE)
    }
    return(deseqds.new)
  }
  dds.fpkm <- dds.filtered
  assay(dds.fpkm, withDimnames = FALSE) <- dds.filtered %>%
    fpkm(robust = FALSE)
  
  dds.tpm <- dds.fpkm
  assay(dds.tpm) <- assay(dds.fpkm) %>%
    sweep(2, colSums(.), '/')
  assay(dds.tpm) <- 1e6 * assay(dds.tpm)
  dds.tpm %<>% perBatchQuantile
  
  # Export
  target_path_counts <- target_folder %>%
    paste("ROSMAP_normalized_log2counts.txt.gz", sep = "")
  dds.tpm %>%
    assay %>%
    write.table(file = gzfile(target_path_counts),
                col.names = NA)
  
  target_path_meta <- target_folder %>%
    paste("ROSMAP_metadata.csv.gz", sep = "")
  dds.tpm %>% colData %>%
    write.table(file = gzfile(target_path_meta),
                col.names = NA)
  
}
# Caret for copy+paste on thinlinc: ^
# Tilde for copy+paste on thinlinc: ~

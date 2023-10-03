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
  jpeg(filename = save_path,
       pointsize = 12,
       height = 600,
       width = 900,
       quality = 95)
  plot(new_plt)
  invisible(dev.off())
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

  #######################
  # CHECK DISTRIBUTIONS #
  #######################
  
  # Plotting function
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
  
  # CPM plot
  plt_path <- paste(plt_folder, "gene_dist.jpeg", sep = "")
  dds.filtered %>%
    checkDistribution(plt_title = "Gene distribution of all samples") %>%
    savePlot(save_path = plt_path) # Plot shows bad overlap and poor normalisation

  # FPKM plot
  fpkm.filtered <- dds.filtered %>%
    fpkm(robust = FALSE)
  plt_path <- paste(plt_folder, "gene_dist_fpkm.jpeg", sep = "")
  log2(1+fpkm.filtered) %>%
    checkDistribution(plt_title = "Gene distribution of all samples after fpkm",
                      xlabel = bquote(log[2](1+FPKM))) %>%
    savePlot(save_path = plt_path)
  
  # TPM plot
  tpm.filtered <- fpkm.filtered %>%
    sweep(2, colSums(.), '/')
  tpm.filtered <- 1e6 * tpm.filtered
  plt_path <- paste(plt_folder, "gene_dist_tpm.jpeg", sep = "")
  log2(1+tpm.filtered) %>%
    checkDistribution(plt_title = "Gene distribution of all samples after tpm",
                      xlabel = bquote(log[2](1+TPM))) %>%
    savePlot(save_path = plt_path)
  
  # Distributions do not have good overlap, use cqn to normalise
  cqn_res <- dds.filtered %>% assay %>%
    cqn(x = gene_data.filtered$percentage_gc_content,
        lengths = gene_data.filtered$gene.length)
  
  checkpoint("cqn_res")
  
  cpm_adj <- (cqn_res$y + cqn_res$offset)
  plt_path <- paste(plt_folder, "gene_dist_cqn.jpeg", sep = "")
  log2(1 + 2^cpm_adj) %>%
    checkDistribution(plt_title = "Gene distribution after cqn") %>%
    savePlot(save_path = plt_path)
  
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
  vsd.filtered.significant_only <- dds.filtered %>%
    getVSDByPadj(min_padj = 0.1)
  vsd.filtered.all_genes <- dds.filtered %>%
    getVSDByPadj(min_padj = 1)
  
  checkpoint("vsd.filtered.all_genes")
  
  # Plot distributions after vst
  plt_path <- paste(plt_folder, "gene_dist_vst_sign.jpeg", sep = "")
  vsd.filtered.significant_only %>% assay %>%
    checkDistribution(plt_title = "Gene distribution of significant genes over all samples after VST",
                      xlabel = bquote(log[2](VST))) %>%
    savePlot(save_path = plt_path)
  
  plt_path <- paste(plt_folder, "gene_dist_vst_all.jpeg", sep = "")
  vsd.filtered.all_genes %>% assay %>%
    checkDistribution(plt_title = "Gene distribution of all samples after VST",
                      xlabel = bquote(log[2](VST))) %>%
    savePlot(save_path = plt_path)
  
  ####################################
  # EXPERIMENTING WITH NORMALISATION #
  ####################################
  
  batches <- vsd.filtered.all_genes$Batch %>% levels
  
  # Sans batch 7
  plt_path <- paste(plt_folder, "gene_dist_vst_all_genes_sans_batch_7.jpeg",
                    sep = "")
  title <- "Gene distribution of all samples after VST sans batch 7"
  sans7_idx <- vsd.filtered.all_genes$Batch != "7"
  vsd.filtered.all_genes[, sans7_idx] %>% assay %>%
    checkDistribution(plt_title = title,
                      xlabel = bquote(log[2](VST))) %>%
    savePlot(save_path = plt_path)
  
  # By batch
  for (batch in batches) {
    
    plt_path <- paste(plt_folder, "gene_dist_vst_all_genes_batch",
                      batch, ".jpeg", sep = "")
    title <- paste("Gene distribution of all samples after VST in batch", batch)
    
    idx <- vsd.filtered.all_genes$Batch == batch
    vsd.filtered.all_genes[, idx] %>% assay %>%
      checkDistribution(plt_title = title,
                        xlabel = bquote(log[2](VST))) %>%
      savePlot(save_path = plt_path)
  }
  
  # CQN on vst sans batch 7
  sans7_idx <- vsd.filtered.all_genes$Batch != "7"
  vst_cqn_res <- 2^(vsd.filtered.all_genes[, sans7_idx] %>% assay) %>%
    cqn(x = gene_data.filtered$percentage_gc_content,
        lengths = gene_data.filtered$gene.length)
  vst_adj <- (vst_cqn_res$y + vst_cqn_res$offset)
  plt_path <- paste(plt_folder, "gene_dist_cqn_post_vst.jpeg", sep = "")
  vst_adj %>%
    checkDistribution(plt_title = "Gene distribution after VST and CQN",
                      xlabel = bquote(log[2](VST))) %>%
    savePlot(save_path = plt_path)
  
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
  
} else {
  #############################################
  # NOT INTERACTIVE, EXPORT NORMALISED COUNTS #
  #############################################
  
  
  # Run DeSeq2
  design(dds.filtered) <- ~ cogdx
  dds.filtered %<>% DESeq(parallel = TRUE)
  
  # Normalise by previous size factors
  vsd.filtered <- dds.filtered %>%
    vst(blind = FALSE)
  
  # Export
  target_path_counts <- target_folder %>%
    paste("ROSMAP_normalized_log2counts.txt.gz", sep = "")
  vsd.filtered %>%
    assay %>%
    write.table(file = gzfile(target_path_counts),
                col.names = NA)
  
  target_path_meta <- target_folder %>%
    paste("ROSMAP_metadata.txt.gz", sep = "")
  dds.supported %>% colData %>%
    write.table(file = gzfile(target_path_meta),
                col.names = NA)
  
}
# Caret for copy+paste on thinlinc: ^
# Tilde for copy+paste on thinlinc: ~

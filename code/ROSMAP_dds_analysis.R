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

working_directory = getwd()

# Needed libraries
library("BiocParallel")
register(MulticoreParam())
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("pheatmap")
library("limma") # removeBatchEffect

# Load the dataset
load(dds_path)

# Pre-filter out very low counts
idx <- dds %>%
  counts %>%
  rowSums(.) >= 10
dds <- dds[idx,]

# Cross-tabulation of cogdx and CERAD
comp_df <- dds %>%
  colData %>%
  as.data.frame
comp_df %<>% .[, c("cogdx", "ceradsc")]
comp_df %<>%
  table %>%
  .[c("1", "2", "4"), ] %>% # exclude the small sample sizes for other cogdx
  prop.table(margin = "cogdx") %>%
  as.data.frame(responseName = "percentcogdx")

comp_df %>% ggplot(aes(x = ceradsc, y = cogdx)) +
  geom_tile(aes(fill = percentcogdx, height = .9)) +
  geom_text(aes(label = percentcogdx %>% round(2)), colour = "white") +
  scale_x_discrete(labels = c("Definite AD", "Probable AD", "Possible AD", "No AD")) +
  scale_y_discrete(labels = c("NCI", "MCI", "AD")) +
  scale_fill_continuous("cogdx%")

# Show that libraryBatch and msex are significant covariates with PCA plot

## Set design as AD or not and perform estimations
design(dds) <- ~ cogdx
dds %<>% DESeq(parallel = TRUE)

## Transform (not with rlog due to speed)
vsd <- dds %>%
  vst(blind = FALSE)

## Show clustering by sex and batch
pcaFromVSD <- function(dst, var = c("libraryBatch", "msex")) {
  pca <- dst %>%
    plotPCA(intgroup = var, returnData = TRUE)
  pca %>% ggplot(aes(PC1, PC2, color = get(var[1]), shape = get(var[2]))) +
    geom_point(size = 3) +
    labs(color = var[1], shape = var[2]) +
    xlab(sprintf("PC1: %2.0f%% variance",100*attr(pca, "percentVar")[1])) +
    ylab(sprintf("PC2: %2.0f%% variance",100*attr(pca, "percentVar")[2]))
}
pcaFromVSD(vsd)

## Show by heatmap the effect of batch and msex (looking at the top x variable genes)

### Show that msex difference is most pronounced among the variable genes
heatmapFromVSD <- function(dst, noGenes, var = c("libraryBatch", "msex")) {
  mat <- head(order(-rowVars(assay(dst))), noGenes) %>%
    assay(dst)[., ]
  mat <- mat - rowMeans(mat)
  
  df <- as.data.frame(colData(dst)[, var])
  pheatmap(mat,
           main = paste("Heatmap of the",noGenes,"most variable genes"),
           annotation_col = df,
           show_colnames = FALSE,
           show_rownames = FALSE)
}
heatmapFromVSD(vsd, 100)

### ...and batch 7 disrupts when the number of genes is higher
heatmapFromVSD(vsd, 2000)

# Because batch 7 has an outsized effect on the data, we remove it
# batch 7 is roughly 8.7% of the data
dds.covar <- dds[, which(dds$libraryBatch != 7)]
dds.covar$libraryBatch %<>% droplevels


#Include msex as covariate in design and estimate new parameters
design(dds.covar) <- ~ msex + cogdx
dds.covar %<>% DESeq(parallel = TRUE)

# Compare new PCA and heatmap
vsd.covar <- dds.covar %>%
  vst(blind = FALSE)
pcaFromVSD(vsd.covar)
heatmapFromVSD(vsd.covar, 1000)

# The heatmap shows that a small number of genes are responsible for msex difference

## Get DE-genes for msex contrast
res.covar.msex <- dds.covar %>%
  results(contrast = c("msex", 0, 1))

## ...and for cogdx contrast
res.covar.cogdx <- dds.covar %>%
  results(contrast = c("cogdx", 4, 1))

## Isolate significant DE set {msex} / {cogdx}
idx_msex_not_AD <- which(
  res.covar.msex$padj < .1 &
  res.covar.cogdx$padj > .1)

## Exclude them from analysis
dds.stripped <- dds.covar[-idx_msex_not_AD, ]

# Redo our analysis for the new dataset
design(dds.stripped) <- ~ cogdx
dds.stripped %<>% DESeq(parallel = TRUE)

# Compare PCA and heatmap again
vsd.stripped <- dds.stripped %>%
  vst(blind = FALSE)

## PCA, of course msex no longer has a significant effect but there is a slight batch effect
pcaFromVSD(vsd.stripped)

## Heatmap, now only RIN and batch seems to be significantly clustered
heatmapFromVSD(vsd.stripped, 1000, c("libraryBatch", "msex", "RIN"))


## We check the PCA with RIN and see a similar effect
pcaFromVSD(vsd.stripped, var = c("RIN", "msex"))

## The RIN <-> batch correlation is reflected in the methods used
### According to the article responsible for RNA-seq (doi:10.1038/s41593-018-0154-9):
### "The libraries were constructed and pooled according to the RIN
### scores such that similar RIN scores would be pooled together.
### Varying RIN scores results in a larger spread of insert sizes
### during library construction and leads to uneven coverage distribution
### throughout the pool."

# Thus we include msex, RIN and batch as covariates and perform one last run-through
dds.final <- dds.covar
design(dds.final) <- ~ msex + libraryBatch + RIN + cogdx
dds.final %<>% DESeq(parallel = TRUE)

# Show the effect of correcting for batch, msex and RIN
vsd.final <- dds.final %>%
  vst(blind = FALSE)

batchMsexRINCorrect <- function(dst, dsgn) {
  mm <- model.matrix(dsgn, data = colData(dst))
  assay(dst) %<>% removeBatchEffect(
    batch = dst$libraryBatch,
    batch2 = dst$msex,
    covariates = dst$RIN,
    design = mm)
  return(dst)
}
vsd.final %<>% batchAndRINCorrect(~ cogdx)

pcaFromVSD(vsd.final, var = c("RIN", "msex"))
pcaFromVSD(vsd.final, var = c("cogdx", "msex"))
pcaFromVSD(vsd.final, var = c("ceradsc_binary", "msex"))
heatmapFromVSD(vsd.final, 2000, var = c("cogdx", "ceradsc_binary"))

# Check results
res.cogdx <- dds.final %>%
  results(contrast = c("cogdx", 4, 1))
res.cogdx %>% summary
res.cogdx %>% DESeq2::plotMA(ylim = c(-1,1),
                       main = "")

## ...for LFC as well (very RAM intensive)
resLFC <- dds.final %>%
  lfcShrink(coef = "cogdx_4_vs_1", type = "apeglm", parallel = TRUE)
resLFC %>% summary
resLFC %>% DESeq2::plotMA(ylim = c(-1,1))
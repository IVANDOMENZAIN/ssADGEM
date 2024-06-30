if (!requireNamespace("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
if (!requireNamespace("viridis", quietly = TRUE)){install.packages("viridis")}
if (!requireNamespace("wacolors", quietly = TRUE)){install.packages("wacolors")}

library("wacolors")
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("viridis")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/plots",sep = "")
dir.create(resultsPath)


counts_matrix <- read_delim(file = 'data/ROSMAP_annotated_samples_counts.txt',delim = '\t', na='NA')
annotation    <- read_delim(file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')

#first, let´s get an overview of the metadata associated to samples
head(annotation)
print(paste("ROSMAP dataset consists of: ", dim(annotation)[[1]],
            " annotated with ", dim(annotation)[[2]]," variables"))

hist(annotation$RIN, main = "RIN distribution")
hist(as.numeric(annotation$Batch), main = "Samples Batch distribution",
     xlab = "Batches")
# Change colors'
colorPallette <- cividis(11)
annotation$Study <- as.factor(annotation$Study)

p <- ggplot(annotation) + theme_bw() +
     geom_bar(aes(x = Study, fill= Study)) + 
     scale_fill_wa_d(wacolors$volcano) 

file_name <- paste(resultsPath,"/ROSMAP_Study.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

p <- ggplot(annotation) + theme_bw() +
  geom_bar(aes(x = msex, fill= msex,xlab("sex"))) + 
  scale_fill_wa_d(wacolors$volcano) 
file_name <- paste(resultsPath,"/ROSMAP_sex.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$educ <- as.numeric(annotation$educ)
p <- ggplot(annotation) + theme_bw() +
     geom_histogram(aes(x = educ)) + 
     scale_fill_wa_d(wacolors$volcano) 
file_name <- paste(resultsPath,"/ROSMAP_educ.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)


annotation$race <- as.factor(annotation$race)
levels(annotation$race) <- c("white","Black","Native american")
p <- ggplot(annotation,aes(x = race)) + theme_bw() +
     geom_bar() + scale_fill_wa_d(wacolors$volcano) 
     file_name <- paste(resultsPath,"/ROSMAP_race.pdf",sep="")
     ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$latinx <- as.factor(annotation$latinx)
levels(annotation$latinx) <- c("Yes","No")
p <- ggplot(annotation,aes(x = latinx)) + theme_bw() +
    geom_bar() + scale_fill_wa_d(wacolors$volcano) 
    file_name <- paste(resultsPath,"/ROSMAP_latinx.pdf",sep="")
    ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$pmi <- as.numeric(annotation$pmi)
p <- ggplot(annotation) + theme_bw() +
     geom_histogram(aes(x = pmi)) + 
     scale_fill_wa_d(wacolors$volcano) 
     file_name <- paste(resultsPath,"/ROSMAP_pmi.pdf",sep="")
     ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

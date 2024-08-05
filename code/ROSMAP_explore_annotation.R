if (!requireNamespace("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
if (!requireNamespace("viridis", quietly = TRUE)){install.packages("viridis")}
if (!requireNamespace("wacolors", quietly = TRUE)){install.packages("wacolors")}

library("wacolors")
library("tidyverse") # Tibble dataframes
library("magrittr") # Piping
library("DESeq2")
library("ggplot2")
library("viridis")
#set wd and create the necessary ones for results
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
base_dir <- getwd()
resultsPath <- paste(base_dir,"/results",sep = "")
dir.create(resultsPath)
resultsPath <- paste(resultsPath,"/plots",sep = "")
dir.create(resultsPath)
#load data
counts_matrix <- read_delim(file = 'data/ROSMAP_annotated_samples_counts.txt',delim = '\t', na='NA')
annotation    <- read_delim(file = 'data/ROSMAP_annotation_samples.txt',delim = '\t', na='NA')
vector 
annotation$AD <- annotation$ceradsc_binary == 'AD' & (annotation$cogdx == 'AD' |  annotation$cogdx == 'AD+') & annotation$braaksc>=5

#first, let´s get an overview of the metadata associated to samples
head(annotation)
print(paste("ROSMAP dataset consists of: ", dim(annotation)[[1]],
            "samples annotated with ", dim(annotation)[[2]]," variables"))

hist(annotation$RIN, main = "RIN distribution")
hist(as.numeric(annotation$Batch), main = "Samples Batch distribution",
     xlab = "Batches")

#understand the distribution of samples regarding different single annotation fields
annotation$Study <- as.factor(annotation$Study)
p <- ggplot(annotation, aes(x = Study, fill= AD)) + theme_bw() +
     geom_bar(position = "stack") + 
     scale_fill_wa_d(wacolors$volcano) 
file_name <- paste(resultsPath,"/ROSMAP_Study.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

names(annotation)[names(annotation) == "msex"] <- "sex"
annotation$sex <- as.factor(annotation$sex)
p <- ggplot(annotation,aes(x = sex, fill= AD)) + theme_bw() +
     geom_bar(position = "stack") + 
     scale_fill_wa_d(wacolors$volcano) 
file_name <- paste(resultsPath,"/ROSMAP_sex.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$educ <- as.numeric(annotation$educ)
p <- ggplot(annotation) + theme_bw() +
     geom_histogram(aes(x = educ,)) + xlab('Years of education') +
     scale_fill_wa_d(wacolors$volcano) 
file_name <- paste(resultsPath,"/ROSMAP_educ.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)


annotation$race <- as.factor(annotation$race)
levels(annotation$race) <- c("white","African american","Native american")
p <- ggplot(annotation,aes(x = race, fill = AD)) + theme_bw() +
     geom_bar(position = "stack") + scale_fill_wa_d(wacolors$volcano) 
     file_name <- paste(resultsPath,"/ROSMAP_race.pdf",sep="")
     ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$latinx <- as.factor(annotation$latinx)
levels(annotation$latinx) <- c("Yes","No")
p <- ggplot(annotation,aes(x = latinx, fill = AD)) + theme_bw() +
    geom_bar(position = "stack") + scale_fill_wa_d(wacolors$volcano) 
    file_name <- paste(resultsPath,"/ROSMAP_latinx.pdf",sep="")
    ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$pmi <- as.numeric(annotation$pmi)
p <- ggplot(annotation) + theme_bw() +
     geom_histogram(aes(x = pmi)) + xlab("post-mortem interval (hours)")
     scale_fill_wa_d(wacolors$volcano) 
     file_name <- paste(resultsPath,"/ROSMAP_pmi.pdf",sep="")
     ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
     
annotation$apoe_genotype <- as.factor(annotation$apoe_genotype)
p <- ggplot(annotation,aes(x = apoe_genotype,fill = AD)) + theme_bw() +
    geom_bar(position = "stack") + scale_fill_wa_d(wacolors$volcano) 
    file_name <- paste(resultsPath,"/ROSMAP_apoe_genotype.pdf",sep="")
    ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
     
annotation$age_death <- as.numeric(annotation$age_death)
p <- ggplot(annotation,aes(x = age_death, fill = AD)) + theme_bw() +
    geom_histogram(position = "stack") + scale_fill_wa_d(wacolors$volcano) 
    file_name <- paste(resultsPath,"/ROSMAP_age_death.pdf",sep="")
    ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
    
annotation$cts_mmse30_lv <- as.numeric(annotation$cts_mmse30_lv)
p <- ggplot(annotation,aes(x = cts_mmse30_lv)) + theme_bw() +
    geom_histogram() + scale_fill_wa_d(wacolors$volcano) +xlab("cognitive test score (last visit)")
    file_name <- paste(resultsPath,"/ROSMAP_cts_mmse30_lv.pdf",sep="")
    ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
    
annotation$braaksc <- as.factor(annotation$braaksc)
p <- ggplot(annotation,aes(x = braaksc, fill = AD)) + theme_bw() +
     geom_bar(position = "stack")+ scale_fill_wa_d(wacolors$volcano) +xlab("Braak stage (histology)")
file_name <- paste(resultsPath,"/ROSMAP_braaksc.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

#ordered_labels <- levels(annotation$ceradsc)
#ordered_labels <- ordered_labels[c(2,3,4,1)]
annotation$ceradsc <- as.factor(annotation$ceradsc)
p <- ggplot(annotation,aes(x = ceradsc, fill = AD)) + theme_bw() +
     geom_bar(position = "stack")+ scale_fill_wa_d(wacolors$volcano) + 
     xlab("CERAD score (plaques)") #+ scale_x_discrete(labels = ordered_labels)
file_name <- paste(resultsPath,"/ROSMAP_ceradsc.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

annotation$ceradsc_binary <- as.factor(annotation$ceradsc_binary)
p <- ggplot(annotation,aes(x = ceradsc_binary, fill = AD)) + theme_bw() +
  geom_bar(position = "stack")+ scale_fill_wa_d(wacolors$volcano) + 
  xlab("CERAD score (binary)") 
file_name <- paste(resultsPath,"/ROSMAP_ceradsc_binary.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

p <- ggplot(annotation,aes(x = cogdx, fill = AD)) + theme_bw() +
  geom_bar(position = "stack")+ scale_fill_wa_d(wacolors$volcano) + 
  xlab("Final cognitive diagnosis") 
file_name <- paste(resultsPath,"/ROSMAP_cogdx.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)

p <- ggplot(annotation,aes(x = dcfdx_lv, fill = AD)) + theme_bw() +
  geom_bar(position = "stack")+ scale_fill_wa_d(wacolors$volcano) + 
  xlab("Final cognitive diagnosis") 
file_name <- paste(resultsPath,"/ROSMAP_cogdx.pdf",sep="")
ggsave(file_name, p, width = 10, height = 10, units = "cm",dpi = 400)
load_data <- function(){
#directory_path %>% setwd()
count_matrix <- "../data/ROSMAP_all_counts_matrix.txt" %>%
  read_delim(delim = "\t") %>%
  column_to_rownames("feature") %>%
  as.matrix()

rnaseq_metadata <- read_csv("../data/ROSMAP_assay_rnaSeq_metadata.csv")
clinical_metadata <- read_csv("../data/ROSMAP_clinical.csv")
biospecimen_metadata <- read_csv("../data/ROSMAP_biospecimen_metadata.csv")
return(list(count_matrix,rnaseq_metadata,clinical_metadata,biospecimen_metadata))
}
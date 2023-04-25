# Descriptive Text


## Filepaths
include filepaths.mk #local to your machine, write your paths in here or change the below variables into paths
human_dir = $(HUMAN_GEM_PATH)
synapse_dir = $(SYNAPSE_DATA_PATH)

## Variables
matlab = matlab -nodisplay -batch

human_files = model/model.mat model/reactions.tsv code/io code/tINIT \
	data/metabolicTasks/metabolicTasks_Essential.txt
human_dep = $(addprefix $(human_dir)/, $(human_files))

ROSMAP_files = ROSMAP_all_counts_matrix.txt ROSMAP_assay_rnaSeq_metadata.csv \
	ROSMAP_clinical.csv ROSMAP_biospecimen_metadata.csv
ROSMAP_dep = $(addprefix $(synapse_dir)/, $(ROSMAP_files))

## Placeholder

build: ;


## Differential Expression Analysis

### Generate Plots
ROSMAP_UMAP.png: data/ROSMAP_dds_processed.rds code/someScript
	touch ROSMAP_UMAP.png

### Test for technical variation
technical_var.png: data/ROSMAP_dds_processed.rds code/someScript
	touch technical_var.png

### Process the data set
ROSMAP_dds_processed.rds: data/ROSMAP_dds.rds code/someScript
	touch ROSMAP_dds_processed.rds

### Create the DESeqDataSet from count matrix and metadata
ROSMAP_dds.rds: $(ROSMAP_dep) code/ROSMAP_fileParser.R
	./code/ROSMAP_fileParser.R $(synapse_dir) data/ROSMAP_dds.rds

## Structural GEM construction

### Construct structural GEMs from a reference model and gene expression data
models.mat: preparedRefModel.mat $(synapse_dir)/ROSMAP_all_counts_matrix.txt code/constructModels.m 
	$(matlab) "code/prepareReferenceModel(\"$(synapse_dir)/ROSMAP_all_counts_matrix.txt\"); exit"

### Prepare the reference model with essential tasks and spontaneous reactions
preparedRefModel.mat: $(human_dep) code/prepareReferenceModel.m
	$(matlab) "code/prepareReferenceModel(\"$(human_dir)\"); exit"


## Functional GEM construction (solving?)

### TODO


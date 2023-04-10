# Descriptive Text


## Filepaths
include filepaths.mk #local to your machine, write your paths in here or change the below variables into paths
human_dir = $(HUMAN_GEM_PATH)
synapse_dir = $(SYNAPSE_DATA_PATH)

## Variables
human_dep = model/model.mat model/reactions.tsv code/io code/tINIT \
	data/metabolicTasks/metabolicTasks_Essential.txt

matlab = matlab -nodisplay -batch

## Placeholder

build: ;


## Differential Expression Analysis

### Generate Plots
ROSMAP_UMAP.png: ROSMAP_dds_processed.rds code/someScript
	touch ROSMAP_UMAP.png

### Test for technical variation
technical_var.png: ROSMAP_dds_processed.rds code/someScript
	touch technical_var.png

### Process the data set
ROSMAP_dds_processed.rds: ROSMAP_dds.rds code/someScript
	touch ROSMAP_dds_processed.rds

### Create the DESeqDataSet from count matrix and metadata
ROSMAP_dds.rds: $(synapse_dir)/[someFiles] code/someScript
	touch ROSMAP_dds.rds

## Structural GEM construction

### Construct structural GEMs from a reference model and gene expression data
models.mat: preparedRefModel.mat $(synapse_dir)/ROSMAP_all_counts_matrix.txt code/constructModels.m 
	$(matlab) "prepareReferenceModel(\"$(synapse_dir)/ROSMAP_all_counts_matrix.txt\"); exit"

### Prepare the reference model with essential tasks and spontaneous reactions
preparedRefModel.mat: $(human_dir)/$(human_dep) code/prepareReferenceModel.m
	$(matlab) "prepareReferenceModel(\"$(human_dir)\"); exit"


## Functional GEM construction (solving?)

### TODO


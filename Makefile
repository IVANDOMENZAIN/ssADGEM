# Descriptive Text


## Filepaths

humanDir = path/to/Human-GEM/ #TODO
synapseDir = /proj/[proj-id]/synapseFiles #TODO
scriptDir = code


## Placeholder

build: ;


## Differential Expression Analysis

### Generate Plots
ROSMAP_UMAP.png: ROSMAP_dds_processed.rds $scriptDir/someScript
	touch ROSMAP_UMAP.png

### Test for technical variation
technical_var.png: ROSMAP_dds_processed.rds $scriptDir/someScript
	touch technical_var.png

### Process the data set
ROSMAP_dds_processed.rds: ROSMAP_dds.rds $scriptDir/someScript
	touch ROSMAP_dds_processed.rds

### Create the DESeqDataSet from count matrix and metadata
ROSMAP_dds.rds: $synapseDir/[someFiles] $scriptDir/someScript
	touch ROSMAP_dds.rds


## Structural GEM construction

### Construct structural GEMs from a reference model and gene expression data
models.mat: preparedRefModel.mat $synapseDir/someFile $scriptDir/someScript 
	touch models.mat

### Prepare the reference model with essential tasks and spontaneous reactions
preparedRefModel.mat: $humanDir/[someFiles] $scriptDir/someScript 
	touch preparedRefModel.mat


## Functional GEM construction (solving?)

### TODO


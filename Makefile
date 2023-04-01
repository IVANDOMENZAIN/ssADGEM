# Descriptive Text

## Filepaths
humanDir = path/to/Human-GEM/ #TODO
rnaDir = path/to/rna #TODO
scriptDir = code

## Placeholder
build: ;

## Differential Expression Analysis
### TODO

## Structural GEM construction

### Construct structural GEMs from a reference model and gene expression data
models.mat: preparedRefModel.mat $transciptionData $GEMscript 
	touch models.mat

### Prepare the reference model with essential tasks and spontaneous reactions
preparedRefModel.mat: $humanDir $prepModelScript
	touch preparedRefModel.mat

## Functional GEM construction (solving?)
### TODO


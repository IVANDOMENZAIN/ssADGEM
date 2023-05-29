% get raven, prepdata, counts, models save location

% add raven

% Load prepared model
prepModel = load(prepModel_path);

% Load the counts table
counts_table = readtable(counts_path);
[~, table_width] = size(counts_table);
preamble_width = 1;
noSamples = table_width - preamble_width;
sampleRange = (1 + preamble_width):table_width;

% Extract information
transcriptStruc.genes = counts_table{:, 1}; % check that these are the gene names
transcriptStruc.tissues = counts_table.Properties.VariableNames(sampleRange); % depends on how readtable works
transcriptStruc.levels = counts_table{:, sampleRange}; % TPM values
transcriptStruc.threshold = 1; % need good number here, what and why???

% Construct the models
models = cell(noSamples, 1);
celltype = [];
hpaData = [];
metabolomicsData = {};
removeGenes = false;
useScoresForTasks = true;

for i = 1:noSamples
    
    models{i} = ftINIT(prepModel, dataStruc.tissues{i}, celltype, ...
        hpaData, transcriptStruc, metabolomicsData, ...
        getHumanGEMINITSteps('1+1'), removeGenes, useScoresForTasks);

end

% Save the models
save(models_savePath, 'models')
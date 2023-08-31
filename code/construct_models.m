function construct_models(prep_model_path, counts_path, hash_salt_path, models_save_dir)

%% Script dependencies

% Find the ssADGEM project dir location
script_path = mfilename('fullpath');
project_path = extractBefore(script_path, "code");

% If no arguments are supplied, assume everything exists in the external or
% data folder
if nargin < 1
    prep_model_path = strjoin([project_path, "data/Human-GEM_reference_model.mat"], "");
    counts_path = strjoin([project_path, "data/ROSMAP_normalized_log2counts.txt"], "");
    hash_salt_path = strjoin([project_path, "data/HASH_salt.txt"], "");
    models_save_dir = strjoin([project_path, "data/"], "");
end

%% Setup

% Load prepared model
prep_model = load(prep_model_path).referenceModel;

% Load the counts table
counts_table = readtable(counts_path, 'PreserveVariableNames', true);
[~, table_width] = size(counts_table);
preamble_width = 1;
no_samples = table_width - preamble_width;
sample_range = (1 + preamble_width):table_width;

% Extract information
transcript_struct.genes = extractBefore(counts_table{:, 1}, "."); % Remove the version number as well
transcript_struct.tissues = counts_table.Properties.VariableNames(sample_range);
transcript_struct.levels = 2.^counts_table{:, sample_range}; % log2 format, otherwise MRN

% Change to CPM?
transcript_struct.levels = 10e6 * transcript_struct.levels ./ sum(transcript_struct.levels,1);

%% Pseudonymisation

% Hash the sample name with salt
% Probably enough for most applications if salt is good enough
salt = fileread(hash_salt_path);

hash_key = table(transpose(transcript_struct.tissues));
hash_key = renamevars(hash_key, "Var1", "Samples");
for i = 1:no_samples

    unhashed_name = [salt transcript_struct.tissues{i}];
    transcript_struct.tissues{i} = string2hash(unhashed_name);

end
hash_key = [hash_key table(transpose(transcript_struct.tissues))];
hash_key = renamevars(hash_key, "Var1",  "Hash");
hash_key_path = strjoin([project_path, "data/HASH_key.csv"], "");
writetable(hash_key, hash_key_path)

%% Construction

% Args
models = cell(no_samples, 1);
celltype = [];
hpaData = [];
metabolomicsData = {}; % TODO: Check if we have this
metsToIgnore = [];
removeGenes = false;
useScoresForTasks = true;

% Start a parallel pool
if isempty(gcp('nocreate'))
    parpool('Processes');
end

% Construct each model
for i = 1:no_samples

    if isfile(strjoin([models_save_dir "ROSMAP_" transcript_struct.tissues{i} ".mat"], ""))
        continue % Don't recreate files already processed
    end

    % We don't want an error in a single model to stop program execution
    try
        model = ftINIT(prep_model, transcript_struct.tissues{i}, celltype, ...
            hpaData, transcript_struct, metabolomicsData, ...
            getINITSteps(metsToIgnore,'1+0'), removeGenes, useScoresForTasks);
    
        save(strjoin([models_save_dir "ROSMAP_" transcript_struct.tissues{i}], ""), "model")
    
    catch
        continue
    end
end

end
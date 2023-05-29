% get inputs: where is the raven, human folder? where should I save 

% add Raven to path
%(do I need to 'check installation' every time?)

% Load the model
humanGEM_path = [humanFolder_path 'model\Human-GEM.mat\'];
human = load("humanGEM_path").model;

% Prep model
addpath([humanFolder_path 'code\tINIT\'])
addpath([humanFolder_path 'code\io\'])
metabolicTasks_path = [humanFolder_path ...
    'data\metabolicTasks\metabolicTasks_Essential.txt'];
reactions_path = [humanFolder_path ...
    'model\reactions.tsv'];
should_convertGenes = false;
prepModel = prepHumanModelForftINIT(human, ...
    should_convertGenes, metabolicTasks_path, reactions_path);

% Save the prepared model for future use
save(prepModel_savePath, 'prepModel')
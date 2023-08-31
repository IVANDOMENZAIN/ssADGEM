function prepare_reference_model(human_dir)

%% Setup

% Find the ssADGEM project dir location
script_path = mfilename('fullpath');
project_path = extractBefore(script_path, "code");
target_path = strjoin([project_path, "data/Human-GEM_reference_model.mat"], "");

% If no directory is supplied, assume symbolic link in 'external' folder
if nargin < 1
    human_dir = strjoin([project_path, "external/human_dir"], "");
end

% Check if default pathdef is writable, if not load the custom pathdef
% which has RAVEN and Gurobi configured
[~, permissions] = fileattrib([matlabroot, '/toolbox/local/pathdef.m']);
if ~permissions.UserWrite

    custom_pathdef_path = strjoin([project_path, ...
        "nobackup/pathdef.m"], "");

    % if RAVEN and Gurobi have not been previously configured, run the
    % install
    if ~isfile(custom_pathdef_path)
        install_path = strjoin([project_path, ...
            "code/install_prereq.m"], "");
        run(install_path)
    end

    addpath( extractBefore(custom_pathdef_path, "pathdef.m") )
    path( pathdef )

end

if isempty(gcp('nocreate'))
    % Start a parallel pool
    parpool('Processes');
end

%% Model prep

% Load the human model
model_path = join([human_dir,"/model/Human-GEM.mat"], "");
human = load(model_path).ihuman;

% Add tINIT and IO functions to path
tINIT_path = join([human_dir,"/code/tINIT"], "");
addpath(tINIT_path);

io_path = join([human_dir,"/code/io"], "");
addpath(io_path);

% Construct the reference model
task_path = join([human_dir, ...
    "/data/metabolicTasks/metabolicTasks_Essential.txt"], "");
reactions_path = join([human_dir, ...
    "/model/reactions.tsv"], "");

should_convertToSymbol = false;

referenceModel = prepHumanModelForftINIT(human, should_convertToSymbol, ...
    task_path, reactions_path);

% Save the reference model
save(target_path, 'referenceModel')

% Remove the added paths again
rmpath(tINIT_path);
rmpath(io_path);

end
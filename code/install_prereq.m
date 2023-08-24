function install_prereq(gurobi_dir, raven_dir)

% Find the ssADGEM project dir location
script_path = mfilename('fullpath');
project_path = extractBefore(script_path, "code");

% If no directory is supplied, assume symbolic link in 'external' folder
if nargin < 1
    gurobi_dir = strjoin([project_path, "external/gurobi_dir"], "");
    raven_dir = strjoin([project_path, "external/raven_dir"], "");
end

try
    % Change dir to gurobi and run gurobi_setup.m
    current_directory = cd(gurobi_dir);
    cd linux64/matlab/
    gurobi_setup
    
    % Change dir to RAVEN and run checkInstall
    cd(raven_dir);
    cd installation
    checkInstallation
    
    % Convince RAVEN to use Gurobi as solver
    setRavenSolver('gurobi')

    % Change back dir
    cd(current_directory)

    % Check if default pathdef is writable, if not save to a custom pathdef
    [~, permissions] = fileattrib([matlabroot, '/toolbox/local/pathdef.m']);
    if ~permissions.UserWrite
        script_path = mfilename('fullpath');
        project_path = extractBefore(script_path, "code");
        custom_pathdef_path = strjoin([project_path, ...
            "nobackup/pathdef.m"], "");
        savepath(custom_pathdef_path)
    % Otherwise save to default
    else
        savepath
    end

catch ME
    % Set user back to working directory in case of error
    cd(current_directory);
    rethrow(ME)
end

end


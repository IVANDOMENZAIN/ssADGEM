function models = models2cell_array(GEMs_dir)
script_path = mfilename('fullpath');
project_path = extractBefore(script_path, "code");
if nargin < 1
    GEMs_dir = strjoin([project_path, "data/ROSMAP_GEMs/"], "");
end
model_files = dir(fullfile(GEMs_dir,'*.mat'));
model_paths = arrayfun(@(x) strjoin([GEMs_dir, x.name], ""), model_files);
models = cell(length(model_paths), 1);
for i = 1:length(model_paths)
    models{i} = load(model_paths(i)).model;
    models{i}.id = model_files(i).name;
end
end
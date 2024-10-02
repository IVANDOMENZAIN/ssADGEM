function sampleIdxs = extractSampleSubset(ADgroup,cluster)
annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');

if nargin<2
    cluster = [];
end
sampleIdxs = annotation.patient(strcmpi(annotation.AD,ADgroup));
%get clusters (GRN) information
if ~isempty(cluster)
    cluster1 = readtable(['../data/samples_GRN_cluster_' num2str(cluster) '.txt'],'Delimiter',',','NumHeaderLines',0);
        %cluster2 = readtable(['../data/samples_GRN_cluster_2.txt'],'Delimiter',',','NumHeaderLines',0);
        %cluster3 = readtable(['../data/samples_GRN_cluster_3.txt'],'Delimiter',',','NumHeaderLines',0);
    cluster1 = strrep(cluster1.Properties.VariableNames,'x','');
    % cluster2 = strrep(cluster2.Properties.VariableNames,'x','');
    % cluster3 = strrep(cluster3.Properties.VariableNames,'x','');
    [~,sampleIdxs] = ismember(cluster1,annotation.patient(strcmpi(annotation.AD,'AD')));
end

%annotation.cluster = zeros(height(annotation),1);
%annotation.cluster(ib(ib>0)) = 1;
%[ia,ib] = ismember(cluster2,annotation.patient(strcmpi(annotation.AD,'AD')));
%annotation.cluster(ib(ib>0)) = 2;
%[ia,ib] = ismember(cluster3,annotation.patient(strcmpi(annotation.AD,'AD')));
%annotation.cluster(ib(ib>0)) = 3;
end

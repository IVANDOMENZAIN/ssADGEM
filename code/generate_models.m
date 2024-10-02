%Get latest version of humanGEM and load xml file with 
git('clone --depth=1 https://github.com/SysBioChalmers/Human-GEM.git')
ihuman = importModel("Human-GEM/model/Human-GEM.xml",false);
mkdir('../results')
%modelVer = ihuman.version;
dateStr  = datetime('today');
geneCounts = readtable('../data/ROSMAP_meanCounts_processed_all.txt','Delimiter','\t');
[~,genes2keep] = ismember(ihuman.genes,table2cell(geneCounts(:,1))); 
geneCounts = geneCounts(genes2keep(genes2keep>0),:);
data_struct = struct();
%filter read by > 1 CPM
ones   = zeros(height(geneCounts),1);
values = table2array(geneCounts(:,2));
values = values/(sum(values)/1E6);
values = log2(values);
idxs   = find(values>0);
%set data structure for ftINIT
%values = table2array(geneCounts(:,2));
values = [values, ones];
data_struct.genes = table2cell(geneCounts(idxs,1));
data_struct.levels = values(idxs,:);
data_struct.tissues = {'ROSMAP' 'ref'};
cd Human-GEM/code/tINIT/
prepDataHumanGEM = prepHumanModelForftINIT(ihuman,false,'../../data/metabolicTasks/metabolicTasks_Essential.txt','../../model/reactions.tsv');
%model1 = ftINIT(prepDataHumanGEM, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
cd ../../../
reducedModel = simplifyModel(prepDataHumanGEM.refModel,false,true,false,true,true);
prepDataHumanGEM.refModel = reducedModel;
[model,~,addedRxnsForTasks, deletedRxnsInINIT] = ftINIT(prepDataHumanGEM, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
model.id = {'Human-GEM pan-ROSMAP'};
save('../models/ROSMAP_GEM.mat','model')
ROSMAP = model;
%generate NCI model
allCounts   = readtable('../data/ROSMAP_counts_processed_all.txt','Delimiter','\t');
annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');
allCounts = allCounts(genes2keep(genes2keep>0),:);
values      = table2array(allCounts(:,1:(end-1)));
[a,b] = size(values);
normVals    = values;
newVals = table2array(geneCounts(:,2));
%values = values/(sum(values)/1E6);
%normalize counts (cpm)
for i=1:b
    normVals(:,i) = values(:,i)./(sum(values(:,i),1)/1E6);
    normVals(normVals(:,i)<1,i) = 0;
end
%normVals = log2(normVals+1);
%normVals = round(normVals,0);
baseModel = ihuman;%prepDataHumanGEM.refModel;
% now build a matrix to compare presence/absence of reactions
NCI     = find(strcmp(annotation.AD,'No_AD'));
NCIcounts = normVals(:,NCI);
logiNCI = logical(NCIcounts);
sumas = sum(logiNCI,2);
idxsSum = find(sumas>=0.5*length(NCI));
tissues = {'NCI' 'ref'};
data_struct.tissues = tissues;
data_struct.genes = geneCounts.gene(idxsSum);
ones = zeros(numel(data_struct.genes),1);
data_struct.levels  = [newVals(idxsSum),ones];
try
    [clusterModel, ~, addedRxnsForTasks, deletedRxnsInINIT]= ftINIT(prepDataHumanGEM,'NCI', [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
    save('../models/ROSMAP_NCI_GEM.mat','clusterModel')
catch
end


%ADpos = find(strcmp(annotation.AD,'AD'));
%cd Human-GEM/code/tINIT/    
patients = allCounts.Properties.VariableNames;
patients = strrep(patients,'x','');
%prepDataHumanGEM2 = prepHumanModelForftINIT(model,false,'../../data/metabolicTasks/metabolicTasks_Essential.txt','../../model/reactions.tsv');
pseudosample = mean(normVals(:,NCI),2);
%[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false, true, false, true, true);
%cd ../../..
compRxnMat= zeros(length(baseModel.rxns), b);
compgnsMat= zeros(length(baseModel.genes), b);
compMetMat= zeros(length(baseModel.mets), b);

compTaskAdd= zeros(length(baseModel.rxns), b);
compINITrmv= zeros(length(baseModel.rxns), b);

%parpool(8)
nonExpressed = cell(1,numel(genes2keep));
for i=1:b
%     %cd Human-GEM/code/tINIT/ 
     patient = annotation.patient{i};
%     disp([num2str(i) ': ' patient])
%     data_struct = struct();
%     data_struct.genes = geneCounts.gene;
%     index   = find(strcmp(patients,patient));
%     ones = zeros(numel(normVals(:,index(1))),1);
%     %data_struct.levels  = [normVals(:,index(1)), normVals(:,NCI)];
%     data_struct.levels  = [normVals(:,index(1)), ones];
% %     nonExpGenes = find(data_struct.levels==0);
% %     nonExpressed{i} = nonExpGenes;
% %     nonExpGenes = data_struct.genes(nonExpGenes);
% %     [~,nonExpGenes] = ismember(nonExpGenes,reducedModel.genes);
% %     sample_model = removeGenes(reducedModel,nonExpGenes,true,true,true);
% %     %cd ../../..
% %     save(['../models/' patient '.mat'],"sample_model")
%     data_struct.tissues = {patient 'NCI'};
%     try
%     [sampleModel, ~, addedRxnsForTasks, deletedRxnsInINIT]= ftINIT(prepDataHumanGEM, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
%     save(['../models/' patient '.mat'],'sampleModel')
    load(['../models_pseudo/' patient '.mat'],'sampleModel')
%     sampleModel.id = {['ROSMAP_' patient]};
    [~,rxnPresence] = ismember(sampleModel.rxns,baseModel.rxns);
    compRxnMat(rxnPresence,i) = 1;
    [~,gnsPresence] = ismember(sampleModel.genes,baseModel.genes);
    compgnsMat(gnsPresence,i) = 1;
    %[~,rxnPresence] = ismember(addedRxnsForTasks,baseModel.rxns);
    %compTaskAdd(rxnPresence,i) = 1;
    %[~,rxnPresence] = ismember(deletedRxnsInINIT,baseModel.rxns);
    %compINITrmv(rxnPresence,i) = 1;
%     catch
%     end
%     clc
end


allCounts   = readtable('../data/ROSMAP_counts_processed_all.txt','Delimiter','\t');
allCounts = allCounts(genes2keep(genes2keep>0),:);
values      = table2array(allCounts(:,1:(end-1)));
[a,b] = size(values);
%values = values/(sum(values)/1E6);
%normalize counts (cpm)
for i=1:b
    normVals(:,i) = values(:,i)./(sum(values(:,i),1)/1E6);
    normVals(normVals(:,i)<1,i) = 0;
end
%normVals = log2(normVals+1);
%normVals = round(normVals,0);
baseModel = ihuman;%pre

%get GRN cluster info
cluster1 = readtable(['../data/samples_GRN_cluster_1.txt'],'Delimiter',',','NumHeaderLines',0);
cluster2 = readtable(['../data/samples_GRN_cluster_2.txt'],'Delimiter',',','NumHeaderLines',0);
cluster3 = readtable(['../data/samples_GRN_cluster_3.txt'],'Delimiter',',','NumHeaderLines',0);

cluster1 = strrep(cluster1.Properties.VariableNames,'x','');
cluster2 = strrep(cluster2.Properties.VariableNames,'x','');
cluster3 = strrep(cluster3.Properties.VariableNames,'x','');

annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');
annotation.cluster = zeros(height(annotation),1);
[ia,ib] = ismember(cluster1,annotation.patient);
annotation.cluster(ib(ib>0)) = 1;
[ia,ib] = ismember(cluster2,annotation.patient);
annotation.cluster(ib(ib>0)) = 2;
[ia,ib] = ismember(cluster3,annotation.patient);
annotation.cluster(ib(ib>0)) = 3;


for i=1:3
    %     %cd Human-GEM/code/tINIT/
    patient_idxs =find(annotation.cluster == i);
    disp(['cluster' num2str(i)])
    data_struct = struct();
    data_struct.tissues = [repelem({'AD'},1,length(patient_idxs)),repelem({'NCI'},1,length(NCI))];
    data_struct.genes = geneCounts.gene;
    ones = zeros(numel(data_struct.genes),1);
    data_struct.levels  = [normVals(:,patient_idxs), normVals(:,NCI)];
    try
        [clusterModel, ~, addedRxnsForTasks, deletedRxnsInINIT]= ftINIT(prepDataHumanGEM,'AD', [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
        save(['../models/AD_GRN_cluster_' num2str(i) '.mat'],'clusterModel')
    catch
    end
end

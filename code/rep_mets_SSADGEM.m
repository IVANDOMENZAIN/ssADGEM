% compRxnMat= zeros(length(baseModel.rxns), b);
% compgnsMat= zeros(length(baseModel.genes), b);
% compTaskAdd= zeros(length(baseModel.rxns), b);
% compINITrmv= zeros(length(baseModel.rxns), b);
 load(['../models_pseudo/ROSMAP_GEM.mat'])
 baseModel = model;
 b = 3;
 repMetsUp = zeros(length(baseModel.mets), b);
 repMetsDn = zeros(length(baseModel.mets), b);
 repMetsAll = zeros(length(baseModel.mets), b);
% nonExpressed = cell(1,numel(genes2keep));
for i=1:3
    data = readtable(['../results/DE_analysis/DE_genes_AD_GRN_cluster_' num2str(i) '.txt'],'Delimiter','\t');
    %load(['../models/AD_GRN_cluster_' num2str(i) '.mat'])
    %model = clusterModel;
    model = baseModel;
    %patient = annotation.patient{i};
    %disp([num2str(i) ': ' patient])
    genes = data.gene;%extractBefore(data.gene,'.');

   %repMets = reporterMetabolites(model,genes,data.padj,true,[],data.log2FoldChange);
   repMets = reporterMetabolites(model,genes,data.padj,true,[],data.log2FC);

   metsUp = repMets(2).mets;
   pVals  = find(repMets(2).metPValues<=0.01);
   if ~isempty(pVals)
   metsUp = metsUp(pVals);
   [~,ib] = ismember(metsUp,baseModel.mets);
   repMetsUp(ib,i) = 1;
   end

   metsDn = repMets(3).mets;
   pVals  = find(repMets(3).metPValues<=0.01);
   if ~isempty(pVals)
   metsDn = metsDn(pVals);
   [~,ib] = ismember(metsDn,baseModel.mets);
   repMetsDn(ib,i) = 1;
   end

    metsAll = repMets(1).mets;
   pVals  = find(repMets(1).metPValues<=0.01);
   if ~isempty(pVals)
   metsAll = metsAll(pVals);
   [~,ib] = ismember(metsAll,baseModel.mets);
   repMetsAll(ib,i) = 1;
   end
   clc
end
mkdir('../results/reporter_metabolites')
%repMetsAll = table(repMetsAll);
writetable(table(repMetsAll),'../results/reporter_metabolites/repMets_GRNclusters_summary_all.txt','Delimiter','\t')
writetable(table(repMetsDn),'../results/reporter_metabolites/repMets_GRNclusters_summary_down.txt','Delimiter','\t')
writetable(table(repMetsUp),'../results/reporter_metabolites/repMets_GRNclusters_summary_up.txt','Delimiter','\t')

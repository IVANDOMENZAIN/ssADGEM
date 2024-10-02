load('Human-GEM/model/Human-GEM.mat')
annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');
ADpos = find(strcmpi(annotation.AD,'AD'));
b = 282;%length(ADpos);
baseModel = ihuman;

compRxnMat= zeros(length(baseModel.rxns), b);
compgnsMat= zeros(length(baseModel.genes), b);
compMetMat= zeros(length(baseModel.mets), b);
for i=1:b
    %data = readtable(['../results/DE_analysis/DE_genes_RM_' num2str(i) '.txt'],'Delimiter','\t');
    patient = annotation.patient{i};
    disp([num2str(i) ': ' patient])
    load(['../models_pseudo/' patient '.mat'],'sampleModel')
    %get model elements presence
    presence = ismember(ihuman.mets,sampleModel.mets);
    compMetMat(presence,i) = 1;
    presence = ismember(ihuman.rxns,sampleModel.rxns);
    compRxnMat(presence,i) = 1;
    presence = ismember(ihuman.genes,sampleModel.genes);
    compgnsMat(presence,i) = 1;
end

rxnPresence = sum(compRxnMat,2);
metPresence = sum(compMetMat,2);
genPresence = sum(compgnsMat,2);
%
figure
set(gca,'FontSize',18)
histogram(rxnPresence(rxnPresence>0),'FaceColor','black')
xlabel('Presence across AD samples')
ylabel('# number of rxns')
figure
set(gca,'FontSize',18)
histogram(metPresence(metPresence>0),'FaceColor','black')
xlabel('Presence across AD samples')
ylabel('# number of mets')
figure
set(gca,'FontSize',18)
histogram(genPresence(genPresence>0),'FaceColor','black')
xlabel('Presence across AD samples')
ylabel('# number of genes')
%
ADpos = find(strcmpi(annotation.AD,'AD'));
ADpos = extractSampleSubset('AD',3);
ADpos = ADpos(ADpos>0);
noADpos = find(~strcmpi(annotation.AD,'AD'));
%Global elements (pan-GEM)
RM_mets = ihuman.mets(metPresence>0);
RM_rxns = ihuman.rxns(rxnPresence>0);
RM_gens = ihuman.genes(genPresence>0);
%Global differences, identify rxns, mets or genes highly differing in
%their presence between AD and non-AD (then also do it by cluster)
AD_rxnPres = sum(compRxnMat(:,ADpos),2);
AD_metPres = sum(compMetMat(:,ADpos),2);
AD_genPres = sum(compgnsMat(:,ADpos),2);
NoAD_rxnPres = sum(compRxnMat(:,noADpos),2);
NoAD_metPres = sum(compMetMat(:,noADpos),2);
NoAD_genPres = sum(compgnsMat(:,noADpos),2);
fc_mets = log2((AD_metPres/length(ADpos))./(NoAD_metPres/length(noADpos)));
fc_rxns = log2((AD_rxnPres/length(ADpos))./(NoAD_rxnPres/length(noADpos)));
fc_gens = log2((AD_genPres/length(ADpos))./(NoAD_genPres/length(noADpos)));

fChange = fc_rxns;
rxnComp = ihuman.rxns;
rxnNames = ihuman.rxnNames;
formulas = constructEquations(ihuman);
subsystems = ihuman.subSystems;
NoADpres  = NoAD_rxnPres/length(noADpos);
ADpres  = AD_rxnPres/length(ADpos);
rxnComp = table(rxnComp,rxnNames,formulas,subsystems,fChange,NoADpres,ADpres);
rxnComp = rxnComp(~isnan(fc_rxns),:);
rxnComp = rxnComp(rxnComp.fChange~=0,:);
rxnComp = sortrows(rxnComp,'fChange','descend');
rxnComp = rxnComp((rxnComp.NoADpres>=0.5),:);
rxnComp = rxnComp(abs(rxnComp.fChange)>0.1,:);

fChange = fc_mets;
metComp = ihuman.mets;
metNames = ihuman.metNames;
components = ihuman.metComps;
NoADpres  = NoAD_metPres/length(noADpos);
ADpres  = AD_metPres/length(ADpos);
metComp = table(metComp,metNames,components,fChange,NoADpres,ADpres);
metComp = metComp(~isnan(fc_mets),:);
metComp = metComp(metComp.fChange~=0,:);
metComp = sortrows(metComp,'fChange','descend');
metComp = metComp((metComp.NoADpres>=0.5),:);
metComp = metComp(abs(metComp.fChange)>0.1,:);

fChange = fc_gens;
genes = ihuman.genes;
shortNAmes = ihuman.geneShortNames;
NoADpres  = NoAD_genPres/length(noADpos);
ADpres  = AD_genPres/length(ADpos);
genComp = table(genes,shortNAmes,fChange,NoADpres,ADpres);
genComp = genComp(~isnan(fc_gens),:);
genComp = genComp(genComp.fChange~=0,:);
genComp = sortrows(genComp,'fChange','descend');
genComp = genComp((genComp.NoADpres>=0.5),:);
genComp = genComp(abs(genComp.fChange)>0.1,:);
%Analyse subsystems
RM_sbsm = ihuman.subSystems(rxnPresence>0);
for i=1:length(RM_sbsm)
    RM_sbsm{i} = char(RM_sbsm{i});
end
unq_RM_sbsm = unique(RM_sbsm);
sbsm_matrix = zeros(numel(unq_RM_sbsm),height(annotation));
tempMat = compRxnMat(find(rxnPresence>0),:);
for i=1:length(unq_RM_sbsm)
    idxs = find(contains(RM_sbsm,unq_RM_sbsm(i)));
    totalRxns = numel(idxs);
    totalSamples = sum(tempMat(idxs,:),1);
    coverage = totalSamples./totalRxns;
    sbsm_matrix(i,:) = coverage;
end
a = mean(sbsm_matrix(:,noADpos),2);
b = mean(sbsm_matrix(:,ADpos),2);
fc = b./a;
log2FC = log2(fc);
mean_cvg_sbsms = table(unq_RM_sbsm,a,b,fc,log2FC);















%cluster = readtable(['../data/samples_GRN_cluster_1.txt'],'Delimiter',',','NumHeaderLines',0);

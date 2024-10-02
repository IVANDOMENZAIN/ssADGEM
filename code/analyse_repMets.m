%analyse_repMets clusters 
load(['../models_pseudo/ROSMAP_GEM.mat'])
all = readtable('../results/reporter_metabolites/repMets_GRNclusters_summary_all.txt','Delimiter','\t');
down = readtable('../results/reporter_metabolites/repMets_GRNclusters_summary_down.txt','Delimiter','\t');
up = readtable('../results/reporter_metabolites/repMets_GRNclusters_summary_up.txt','Delimiter','\t');

matrix = table2array(all);
sumas = sum((matrix),2);
x = find(sumas==1);
mets = model.mets(x);
names = model.metNames(x);
compartments = model.compNames(model.metComps(x));
overlapTable = table(mets,names,compartments);
c1 = model.mets(find(matrix(:,1)==1));
c2 = model.mets(find(matrix(:,2)==1));
c3 = model.mets(find(matrix(:,3)==1));

cc1 = model.metNames(find(matrix(:,1)==1));
cc2 = model.metNames(find(matrix(:,2)==1));
cc3 = model.metNames(find(matrix(:,3)==1));

ccm1 = model.compNames(model.metComps(find(matrix(:,1)==1)));
ccm2 = model.compNames(model.metComps(find(matrix(:,2)==1)));
ccm3 = model.compNames(model.metComps(find(matrix(:,3)==1)));
newTable = table(cc3,ccm3);
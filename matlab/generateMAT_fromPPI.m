
clear all;
%what file has all the values
labledCSV='../PPI/data.csv';

%%%%%%%load features
%from csv
genes.feats  = csvread(labledCSV,1,1);
%set lables
genes.labels = genes.feats(:,end);
%remove lables from feats
genes.feats  = genes.feats(:,1:end-1);


save('PPI.mat', 'genes');

fprintf('+ %i\t - %i (%i)\n', length(find(genes.labels==1)), length(find(genes.labels==0)),length(genes.labels));



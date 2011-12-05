clear all;

load 'genes.mat';
%genes = 
%  feats: [5418x11 double]
%  labels: [5418x1 double]
%  featnames: {1x11 cell}
%  names: {1x5418 cell}

%c         = cvpartition(genes.labels,'holdout',.3);
%unlabeled = find(c.training);
%heldout   = find(c.test);
%labeled   = [];


%%%
%random initial (10%)
%%
avalableIdx = [1:length(genes.feats)];
trainIdx    = randsample( avalableIdx, floor(length(avalableIdx)/10) );
avalableIdx(trainIdx) = [];


svm = svmtrain( genes.feats(trainIdx,:), genes.labels(trainIdx,:) );

results = svmclassify(svm,genes.feats);

%stats
TP = length( find( results( find(genes.labels == 1) ) ==1) );
TN = length( find( results( find(genes.labels == 0) ) ==0) );

FP = length( find( results( find(genes.labels == 0) ) ==1) );
FN = length( find( results( find(genes.labels == 1) ) ==0) );


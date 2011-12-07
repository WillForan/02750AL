function initial
    %clear all;
    global genes;
    global avalableIdx;
    global trainIdx;
    global correct;

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

    %uses trainIdx and avalableIdx
    i=1;
    randtrain(i)

    %%%% random
    batch=100;
    while(avalableIdx(1))
	remain      = length(avalableIdx);
	if(remain < batch); batch=remain; end; 

	removeIdxIdx = randsample( length(avalableIdx), batch );
	trainIdx=[trainIdx avalableIdx(removeIdxIdx)];
	avalableIdx(removeIdxIdx) = [];


	%length(unique(trainIdx)) == length(trainIdx)


	i=i+1;
	randtrain(i);
    end
end

function randtrain(i)
    global avalableIdx; global trainIdx;
    global genes;
    global correct;
    %train on random
    svm = svmtrain( genes.feats(trainIdx,:), genes.labels(trainIdx,:) );

    %get total results
    results = svmclassify(svm,genes.feats);

    %stats
    TP = length( find( results( find(genes.labels == 1) ) ==1) );
    TN = length( find( results( find(genes.labels == 0) ) ==0) );
    FP = length( find( results( find(genes.labels == 0) ) ==1) );
    FN = length( find( results( find(genes.labels == 1) ) ==0) );
    correct(i)=(TP+TN)/length(genes.labels)*100;
    fprintf('[*] %i - %i (%i)\n', i, length(trainIdx), length(unique(trainIdx)));
    fprintf('\tTP\tTN\t\tFP\tFN\tcorrect\n\t%i\t%i\t\t%i\t%i\t%2.1f\n\n', TP,TN,FP,FN,correct(i));
end

function defineGlobals
    global genes;
    global avalableIdx;
    global trainIdx;
    global correct;
end

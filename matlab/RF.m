function initial
    %clear all;
    globals;

    load 'genes.mat';
    validationFolds=3;
    path(path,'randomforest-matlab-read-only/RF_Class_C/');

    %genes = 
    %  feats: [5418x11 double]
    %  labels: [5418x1 double]
    %  featnames: {1x11 cell}
    %  names: {1x5418 cell}

    %c         = cvpartition(genes.labels,'holdout',.3);
    %unlabeled = find(c.training);
    %heldout   = find(c.test);
    %labeled   = [];
    c = cvpartition(genes.labels,'k',validationFolds);
    for f=1:c.NumTestSets
	RF               = classRF_train(genes.feats(c.training(f),:),genes.labels(c.training(f)),500);
	[results, votes] = classRF_predict(genes.feats(c.test(f)),RF);

	%positive index of labels in test
        posIdx = find(genes.labels(c.test(f)) == 1);
        negIdx = find(genes.labels(c.test(f)) == 0);

	TP = length(  find( results(posIdx) ==1 )  );
	TN = length(  find( results(negIdx) ==0 )  );
	FP = length(  find( results(negIdx) ==1 )  );
	FN = length(  find( results(posIdx) ==0 )  );

	accuracy    = (TP+TN)/length(c.test(f));
	sensitivity =      TP/(TP+FN);
	specificity =      TP/(TN+FP);
	FPR	    =      FP/(FP+TN);

	fprintf(['*\tTP\tTN\t\tFP\tFN\tacc\tsens\tspec\tFPR\n',     ...
	         '\t%i\t%i\t\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\n\n'], ...
		TP,TN,FP,FN,accuracy,sensitivity, specificity,FPR );

	%label	predict	diff	|     votes
	% 1	1	0   TP	| >0  more for 1   P
	% 0	0	0   TN	| <0  more for 0   N
	% 1	0	1   FN	|
	% 0	1	-1  FP	|
	votediff = [genes.labels(c.test(f)) - results, votes(:,1) - votes(:,2)];
	%the avg of how many more votes went positive for correctly classified
	FNdiff   = sum( votediff(find(votediff(:,1)==1),2) )/FN;
	FPdiff   = sum( votediff(find(votediff(:,1)==-1),2) )/FP;
	%first is 0 (no diff==match), second is which side of 0
	TPdiff   = sum( votediff(find( votediff(find(votediff(:,1)==0),2) > 0 ),2) )/TP;
	TNdiff   = sum( votediff(find( votediff(find(votediff(:,1)==0),2) < 0 ),2) )/TN;

	fprintf(['avg\ttpdif\ttndiff\t\tfpdiff\tfndiff\n', ...
	         '\t%2.1f\t%2.1f\t\t%2.1f\t%2.1f\n\n\n'],  ...
		TPdiff,TNdiff,FPdiff,FNdiff);



    end

    %err = zeros(CVO.NumTestSets,1);
    %for i = 1:CVO.NumTestSets
    %    trIdx = CVO.training(i);
    %    teIdx = CVO.test(i);
    %    ytest = classify(meas(teIdx,:),meas(trIdx,:),...
    %    species(trIdx,:));
    %    err(i) = sum(~strcmp(ytest,species(teIdx)));
    %end
    %cvErr = sum(err)/sum(CVO.TestSize)
    return

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
    globals;
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


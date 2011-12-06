function initial
    %clear all;
    globals;

    load 'genes.mat';
    validationFolds=5;
    numTrees=50;
    path(path,'randomforest-matlab-read-only/RF_Class_C/');

    %genes = 
    %  feats: [5418x11 double]
    %  labels: [5418x1 double]
    %  featnames: {1x11 cell}
    %  names: {1x5418 cell}

    c = cvpartition(genes.labels,'k',validationFolds);
    r=struct;

    %for each fold
    for f=1:c.NumTestSets
	RF               = classRF_train(genes.feats(c.training(f),:),genes.labels(c.training(f)),numTrees);
	[results, votes] = classRF_predict(genes.feats(c.test(f)),RF);

	%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%

	%positive index of labels in test
        posIdx = find(genes.labels(c.test(f)) == 1);
        negIdx = find(genes.labels(c.test(f)) == 0);

	r.TP(f) = length(  find( results(posIdx) ==1 )  );
	r.TN(f) = length(  find( results(negIdx) ==0 )  );
	r.FP(f) = length(  find( results(negIdx) ==1 )  );
	r.FN(f) = length(  find( results(posIdx) ==0 )  );

	r.accuracy(f)	 = (r.TP(f)+r.TN(f))/c.TestSize(f);
	r.sensitivity(f) =           r.TP(f)/( r.TP(f) + r.FN(f) );
	r.specificity(f) =           r.TP(f)/( r.TN(f) + r.FP(f) );
	r.FPR(f)	 =           r.FP(f)/( r.FP(f) + r.TN(f) );


	%label	predict	diff	|     votes
	%			| class(1) - class(0)
	% 1	1	0   TP	| >0  more for 1   P
	% 0	0	0   TN	| <0  more for 0   N
	% 1	0	1   FN	|
	% 0	1	-1  FP	|
	votediff = [genes.labels(c.test(f)) - results, votes(:,2) - votes(:,1)];
	%the avg of how many more votes went positive for correctly classified
	r.FNdiff(f)   = sum( votediff(find(votediff(:,1)==1),2) )/r.FN(f);
	r.FPdiff(f)   = sum( votediff(find(votediff(:,1)==-1),2) )/r.FP(f);
	%trues are (1) == 0 (no diff==match)
	%find which side of 0 trutth si on to deterim if TN or TP
	trues         = votediff(find(votediff(:,1)==0),:);
	r.TPdiff(f)   = sum( trues( find( trues(:,2) > 0 ),2)) / r.TP(f);
	r.TNdiff(f)   = sum( trues( find( trues(:,2) < 0 ),2)) / r.TN(f);

    end

    names=fieldnames(r);
    for n=1:length(names)
	fprintf('\t%s\t%2.2f\n', ...
	    names{n},           ...
	    sum(getfield(r,names{n}))/c.NumTestSets);
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


    %%%%%%%%%
    %for tree
    %%%%%%%%%
    votediff = [genes.labels(c.test(f)) - results, votes(:,1) - votes(:,2)];
    %the avg of how many more votes went positive for correctly classified
    FNdiff  = sum( votediff(find(votediff(:,1)==1),2) )/FN;
    FPdiff  = sum( votediff(find(votediff(:,1)==-1),2) )/FP;
    %first is 0 (no diff==match), second is which side of 0
    TPdiff(f)   = sum( votediff(find( votediff(find(votediff(:,1)==0),2) > 0 ),2) )/TP;
    TNdiff(f)   = sum( votediff(find( votediff(find(votediff(:,1)==0),2) < 0 ),2) )/TN;
    %%%%%%%%%

    correct(i)=(TP+TN)/length(genes.labels)*100;
    fprintf('[*] %i - %i (%i)\n', i, length(trainIdx), length(unique(trainIdx)));
    fprintf('\tTP\tTN\t\tFP\tFN\tcorrect\n\t%i\t%i\t\t%i\t%i\t%2.1f\n\n', TP,TN,FP,FN,correct(i));
end


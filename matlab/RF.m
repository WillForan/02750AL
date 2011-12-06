function initial
    %clear all;
    globals;

    load 'genes.mat';
    validationFolds=5;
    numTrees=50;
    BatchSize=100;
    path(path,'randomforest-matlab-read-only/RF_Class_C/');

    %genes = 
    %  feats: [5418x11 double]
    %  labels: [5418x1 double]
    %  featnames: {1x11 cell}
    %  names: {1x5418 cell}

    c = cvpartition(genes.labels,'k',validationFolds);
    r=struct;

    %accuracy for each fold and each interval step
    maxIterations  = ceil(max(c.TrainSize)/BatchSize);
    stats.accuracy = zeros(c.NumTestSets, maxIterations );
    numInstStart   = floor(min(c.TrainSize)/10); %always start with 10% of the smallest training set

    %for each fold
    for f=1:c.NumTestSets
	
	fprintf('[*][*] %i/%i\n',f,c.NumTestSets);

	%%%%%%%%%%%%%%%%%%
	% ACTIVE LEARNING
	%%%%%%%%%%%%%%%%%%

	%%%
	%random initial (10%)
	%%
	avalableIdx  = find(c.training(f)==1)';			%all in the training set are avalable to be introduced to AL traning
	removeIdxIdx = randsample( length(avalableIdx), ...
			        numInstStart);			%pick 10% to remove (move from avalable to AL train)
	trainIdx     = avalableIdx(removeIdxIdx); 		%add the indexes to trainIdx
	avalableIdx(removeIdxIdx)=[]; 				%remove avalability

	%%%%%%
	% start AL iterations
	%%%%%%
	i=1;
	trainRF(f,i); 

	%%%% contine until there are no points left
	while(length(avalableIdx)>0)
	    i=i+1;
	    updateTrainIdxRND;
	    trainRF(f,i);
	end
    end


    TP = sum(stats.random.TP ) ./ c.NumTestSets;
    TN = sum(stats.random.TN ) ./ c.NumTestSets;
    FP = sum(stats.random.FP ) ./ c.NumTestSets;
    FN = sum(stats.random.FN ) ./ c.NumTestSets;

    fig=figure;
    %plot ([1:maxIterations] .* BatchSize,sum(stats.random.accuracy) ./ c.NumTestSets);
    plot (...
          [ numInstStart:BatchSize:BatchSize*length(stats.random.TP)+BatchSize], ...
          [(TP+TN) ./ (TP+TN+FP+FN); TP./(TP+FN); TP./(TN+FP); FP./(FP+TN)  ]...
	 );
    legend({'Accuracy','Sensitivity','Specificity','FDR'});
    title(['Random Forest - Random selection - Batch ' num2str(BatchSize)]);
    xlabel('Instances Seen'); ylabel(''); 
    hgexport(fig,'img/random')


    %names=fieldnames(r);
    %for n=1:length(names)
    %    fprintf('\t%s\t%2.2f\n', ...
    %        names{n},            ...
    %        sum(getfield(r,names{n}))/c.NumTestSets);
    %end


end

function trainRF(f,i)
    globals;

    %train on selected features
    clear model;
    model = classRF_train( genes.feats(trainIdx,:), genes.labels(trainIdx,:),numTrees );

    %get total results
    results=zeros(1,c.TestSize(f));
    votes=zeros(2,c.TestSize(f));
    [results,votes] = classRF_predict(genes.feats(c.test(f)),model);


    %%%%%%%%%
    %for tree
    %%%%%%%%%

    posIdx = find(genes.labels(c.test(f)) == 1);
    negIdx = find(genes.labels(c.test(f)) == 0);
    stats.random.TP(f,i) = length(  find( results(posIdx) ==1 )  );
    stats.random.TN(f,i) = length(  find( results(negIdx) ==0 )  );
    stats.random.FP(f,i) = length(  find( results(negIdx) ==1 )  );
    stats.random.FN(f,i) = length(  find( results(posIdx) ==0 )  );

    fprintf('[*] %i - %i\n', i, length(trainIdx));
    %votediff = [genes.labels(c.test(f)) - results, votes(:,1) - votes(:,2)];
    %trues         = votediff(find(votediff(:,1)==0),:);
    %stats.random.accuracy(f,i)=length(trues)/c.TestSize(f);
    %fprintf('[*] %i - %i\t%.3f\n', i, length(trainIdx),stats.random.accuracy(f,i));
end

function updateTrainIdxRND
    globals;
    remain = length(avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 

    removeIdxIdx = randsample( length(avalableIdx), batch );
    trainIdx=[trainIdx avalableIdx(removeIdxIdx)];
    avalableIdx(removeIdxIdx) = [];
end

function updateTrainIdxVOTES
    globals;
    %set batch
    remain = length(avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 

    removeIdxIdx = randsample( length(avalableIdx), batch );
    trainIdx=[trainIdx avalableIdx(removeIdxIdx)];
    avalableIdx(removeIdxIdx) = [];
end

function treeJunk
	%%%%%%%%%%%%%%%%%%
	RF               = classRF_train(genes.feats(c.training(f),:),genes.labels(c.training(f)),numTrees);
	[results, votes] = classRF_predict(genes.feats(c.test(f)),RF);

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

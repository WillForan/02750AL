function initial
    %clear all;
    globals;

    load 'genes.mat';
    validationFolds=3;
    numTrees=30;
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

        %count the number of differentially expressed genes
	%allPosi(f) = length(find(genes.labels(avalableIdx)==1));

	%%%%%%
	% start AL iterations
	%%%%%%
	i=1;

	%%%% contine until there are no points left
	while(length(avalableIdx)>0)
	    trainRF(f,i);
	    classifyRF_ontest(f,i);
	    classifyRF_ontrain(f,i);
	    i=i+1;
	    updateTrainIdxMostVotes
	    %updateTrainIdxRND;
	    %updateTrainIdxMostVotes
	end
    end


    aTP = sum(stats.random.ontrain.TP ) ./ c.NumTestSets;
    aTN = sum(stats.random.ontrain.TN ) ./ c.NumTestSets;
    aFP = sum(stats.random.ontrain.FP ) ./ c.NumTestSets;
    aFN = sum(stats.random.ontrain.FN ) ./ c.NumTestSets;

    for f=1:c.NumTestSets; TPr(f,:) = stats.random.ontrain.TP(f,:) ./ length(find(genes.labels(c.training(f))==1)); end;
    TPr = sum(TPr) ./ c.NumTestSets;

    bTP = sum(stats.random.ontest.TP ) ./ c.NumTestSets;
    bTN = sum(stats.random.ontest.TN ) ./ c.NumTestSets;
    bFP = sum(stats.random.ontest.FP ) ./ c.NumTestSets;
    bFN = sum(stats.random.ontest.FN ) ./ c.NumTestSets;

    fig=figure;
    %plot ([1:maxIterations] .* BatchSize,sum(stats.random.accuracy) ./ c.NumTestSets);
    plot (...
          [ numInstStart:BatchSize:BatchSize*(length(stats.random.ontest.TP)-1)+numInstStart], ...
          [ TPr
	    (aTP + aTN)./(aTP + aTN + aFP + aFN)
	    (bTP + bTN)./(bTP + bTN + bFP + bFN) ...
	    ] ...
         );

	    %TP./(TP+FN) ... Sensitivity
	    %TP./(TN+FP) ... Sepecificity
	    %FP./(FP+TN) ... 
    legend({'TruePos:Total','Train Accuracy','Test Accuracy'},'Location', 'NorthOutside');
    title(['Random Forest - Most Votes selection - Batch Size ' num2str(BatchSize) '- ' num2str(validationFolds) ' Folds']);
    xlabel('Instances Seen'); ylabel(''); 
    hgexport(fig,'img/RF-MV')


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
    model =[];
    model = classRF_train( genes.feats(trainIdx,:), genes.labels(trainIdx,:),numTrees );
end

function classifyRF_ontest(f,i)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    results = classRF_predict(genes.feats(c.test(f)),model);


    %%%%%%%%%
    %for tree
    %%%%%%%%%

    posIdx = find(genes.labels(c.test(f)) == 1);
    negIdx = find(genes.labels(c.test(f)) == 0);
    stats.random.ontest.TP(f,i) = length(  find( results(posIdx) ==1 )  );
    stats.random.ontest.TN(f,i) = length(  find( results(negIdx) ==0 )  );
    stats.random.ontest.FP(f,i) = length(  find( results(negIdx) ==1 )  );
    stats.random.ontest.FN(f,i) = length(  find( results(posIdx) ==0 )  );

    fprintf('[** test] %i - %i\n', i, length(trainIdx));
    %votediff = [genes.labels(c.test(f)) - results, votes(:,1) - votes(:,2)];
    %trues         = votediff(find(votediff(:,1)==0),:);
    %stats.random.accuracy(f,i)=length(trues)/c.TestSize(f);
    %fprintf('[*] %i - %i\t%.3f\n', i, length(trainIdx),stats.random.accuracy(f,i));
end
function classifyRF_ontrain(f,i)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    stats.random.ontrain.votes=[];%zeros(1,genes.feats(avalableIdx));
    [results,votes] = classRF_predict(genes.feats(avalableIdx,:),model);


    %
    % Idx of avalableIdx ranked by vote difference
    [~,stats.random.ontrain.votes]=sort(votes(:,1) - votes(:,2));



    posIdx = find(genes.labels(avalableIdx) == 1);
    negIdx = find(genes.labels(avalableIdx) == 0);

    stats.random.ontrain.TP(f,i) = length(  find( results(posIdx) ==1 )  );
    stats.random.ontrain.TN(f,i) = length(  find( results(negIdx) ==0 )  );
    stats.random.ontrain.FP(f,i) = length(  find( results(negIdx) ==1 )  );
    stats.random.ontrain.FN(f,i) = length(  find( results(posIdx) ==0 )  );

    fprintf('[** train] %i - %i\n', i, length(trainIdx));
    return
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

function updateTrainIdxMostVotes
    globals;
    %set batch
    remain = length(avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 

    removeIdxIdx = stats.random.ontrain.votes(end-batch:end);
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

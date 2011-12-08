function RF
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
    %load('partion.mat'); %loads c
    r=struct;

    %%%%%%%%%%%%%%%%%%%
    %define methods
    %%%%%%%%%%%%%%%%%%%
    loadMethods; %e.g. AL.RFLV AL.RFRandom
    usemethods={'RFRandom'};
    %usemethods={'RFMostVotes',  'RFLV'};
    %usemethods=fieldnames(AL)';

    %accuracy for each fold and each interval step
    maxIterations  = ceil(max(c.TrainSize)/BatchSize);
    %AL.accuracy = zeros(c.NumTestSets, maxIterations );
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


	for alType=usemethods
	    alType=alType{1};
	    disp(alType);
	    
	    AL.(alType).ontrain.avalableIdx=avalableIdx;
	    AL.(alType).ontrain.trainIdx=trainIdx;

	    %count the number of differentially expressed genes
	    %allPosi(f) = length(find(genes.labels(avalableIdx)==1));

	    %%%%%%
	    % start AL iterations
	    %%%%%%
	    i=1;

	    %%%% contine until there are no points left
	    while(length(AL.(alType).ontrain.avalableIdx)>0)

		%%%%%MODEL
		%populate model (global) with trainIdx
		AL.(alType).trainer(f,i,AL.(alType).ontrain.trainIdx);
		%classify this model on what is left out for testing
		AL.(alType).ontest  = AL.(alType).testCly(f,i, AL.(alType).ontest);
		%classify this model on what we've seen so far
		AL.(alType).ontrain = AL.(alType).trainCly(f,i,AL.(alType).ontrain);


		%%%%Instance selection for next training round
		%we're on the next step now
		i=i+1;
		%pick what to look at next
		AL.(alType).ontrain = AL.(alType).updater(AL.(alType).ontrain);
		%%choose from
		%%updateTrainIdxRND;
		%%updateTrainIdxMostVotes
	    end
	end
    end

    for alType=usemethods
	%add values used, and save
	alType=alType{1};
	AL.(alType).validationFolds=validationFolds;
	AL.(alType).BatchSize=BatchSize;
	AL.(alType).numTrees=numTrees;
	AL.(alType).training=arrayfun(@(x)(c.training(x)), [1:validationFolds],'UniformOutput',0);
	AL.(alType).numInstStart=numInstStart;
	save(['alResults/' (alType) '.mat'],'-struct', 'AL', alType);
    end
    %save partion stuff
    save('partion.mat','c');





end

function trainRF(f,i,trainIdxs)
    clear model
    globals;

    %train on selected features
    model = classRF_train( genes.feats(trainIdxs,:), genes.labels(trainIdxs,:),numTrees );
end

function ontest=classifyRF_ontest(f,i,ontest)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    results = classRF_predict(genes.feats(c.test(f)),model);


    %%%%%%%%%
    %%for tree
    %%%%%%%%%

    posIdx = find(genes.labels(c.test(f)) == 1);
    negIdx = find(genes.labels(c.test(f)) == 0);

    ontest.TP(f,i) = length(  find( results(posIdx) ==1 )  );
    ontest.TN(f,i) = length(  find( results(negIdx) ==0 )  );
    ontest.FP(f,i) = length(  find( results(negIdx) ==1 )  );
    ontest.FN(f,i) = length(  find( results(posIdx) ==0 )  );

    fprintf('[** test] %i\n', i);
    %votediff = [genes.labels(c.test(f)) - results, votes(:,1) - votes(:,2)];
    %trues         = votediff(find(votediff(:,1)==0),:);
    %AL.RFMostVotes.accuracy(f,i)=length(trues)/c.TestSize(f);
    %fprintf('[*] %i - %i\t%.3f\n', i, length(trainIdx),AL.RFMostVotes.accuracy(f,i));
end
function ontrain = classifyRF_ontrain(f,i,ontrain)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    Idx=[ontrain.trainIdx ontrain.avalableIdx];
    results = classRF_predict(genes.feats(Idx,:),model);



    posIdx = find(genes.labels(Idx) == 1);
    negIdx = find(genes.labels(Idx) == 0);

    ontrain.TP(f,i) = length(  find( results(posIdx) ==1 )  );
    ontrain.TN(f,i) = length(  find( results(negIdx) ==0 )  );
    ontrain.FP(f,i) = length(  find( results(negIdx) ==1 )  );
    ontrain.FN(f,i) = length(  find( results(posIdx) ==0 )  );

    fprintf('[** train] %i - %i\n', i, length(ontrain.trainIdx));
    return
end

function ontrain = updateTrainIdxRND(ontrain)
    globals;
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 

    removeIdxIdx = randsample( length(ontrain.avalableIdx), batch );
    ontrain.trainIdx=[ontrain.trainIdx ontrain.avalableIdx(removeIdxIdx)];
    ontrain.avalableIdx(removeIdxIdx) = [];
end

function ontrain = updateTrainIdxMostVotes(ontrain)
    globals;

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote difference
    [results,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    [~,votes]       = sort(abs(votes(:,1) - votes(:,2))); %the bigger the difference, the better

    %set batch
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    removeIdxIdx     = votes(end-batch:end);
    ontrain.trainIdx = [ontrain.trainIdx ontrain.avalableIdx(removeIdxIdx)];
    ontrain.avalableIdx(removeIdxIdx) = [];
end
function ontrain = updateTrainIdxMostPositiveVotes(ontrain)
    globals;

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote difference
    [results,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    [~,votes]       = sort(votes(:,1) - votes(:,2)); %positive votes are true class

    %set batch
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    removeIdxIdx     = votes(end-batch:end);
    ontrain.trainIdx = [ontrain.trainIdx ontrain.avalableIdx(removeIdxIdx)];
    ontrain.avalableIdx(removeIdxIdx) = [];
end
function ontrain = updateTrainIdxLeastVotes(ontrain)
    globals;

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote difference
    [results,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    [~,votes]       = sort(abs(votes(:,1) - votes(:,2)),'descend'); %sort smallest votedifference first

    %set batch
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    removeIdxIdx     = votes(end-batch:end);
    ontrain.trainIdx = [ontrain.trainIdx ontrain.avalableIdx(removeIdxIdx)];
    ontrain.avalableIdx(removeIdxIdx) = [];
end


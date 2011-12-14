% for i=1:5; for d={'PPI' 'PF'}; runAL({'RFKMeans' 'RFConfusion'  'RFMostVotes' 'RFMostPositiveVotes' 'SMORND' 'RFRandom' 'RFEntropy'},d{1},i); end; end;     
function runAL(usemethods, whichdata, whichiteration)
    %clear all;
    globals;

    validationFolds=5;
    numTrees=30;
    numInstances=15; %determines batchSize

    %Default to genes
    %load 'genes.mat';ResultsFolder='alResults';
    %give genes.feats and genes.labels

    %Can also load PPI
    %load 'PPI.mat'; ResultsFolder='PPIResults';
    % still has genes.feats and genes.labels structure

    if(exist('whichdata'))
	switch whichdata
	case {'PPI' 'ppi'}
	    loadPPI;
	case {'oldPF','oldgenes','PFnogc'}
	    loadGenes_nogc;
	case {'PF','genes'}
	    loadGenes;
	otherwise
	    loadGenes;
	end
    else
	loadGenes;
    end

    if(~exist('whichiteration'))
	whichiteration='';
    else
	whichiteration=num2str(whichiteration);
    end


    %%%%% TESTING
    %ResultsFolder='Tests'
    %validationFolds=2;
    %numTrees=15;
    %numInstances=5; %determines batchSize
    %%%%%

    %load Random Forest MEX (unstable!but svn@46 doesn't compile)
    path(path,'randomforest-matlab-read-only/RF_Class_C/');


    %parition data for kfold cross-validation
    c = cvpartition(genes.labels,'k',validationFolds);

    %BatchSize is given by the numInstances we want and the size of what we have to work with
    BatchSize     = floor(min(c.TrainSize)/numInstances); %floor to sqeeze the last incomplete batch in
    maxIterations = ceil(max(c.TrainSize)/BatchSize);

    %%dont do this, we want different (random) partions
    %load('partion.mat'); %loads c 
    %%could be useful to keep trails compareable

    %results structure -- not used anymore
    %r=struct;

    %%%%%%%%%%%%%%%%%%%
    %define Active Learning methods
    %%%%%%%%%%%%%%%%%%%
    loadMethods; %e.g. AL.RFLV AL.RFRandom
    
    if(~exist('usemethods'))
	%usemethods={'RFRandom'};
	%usemethods={'RFMostVotes',  'RFLV'};
	usemethods=fieldnames(AL)';
    end

    %accuracy for each fold and each interval step
    %AL.accuracy = zeros(c.NumTestSets, maxIterations );

    %% if starting with 10% random, should use this -- not used
    %numInstStart   = floor(min(c.TrainSize)/10); %always start with 10% of the smallest training set

    %%%%%%%
    %for each fold
    %   for each AL type
    %       until all instances are classified
    %%%%%%%
    for f=1:c.NumTestSets
	fold=f;%used f inside functions then decided easier for fold to be global. oops
	
	fprintf('[*][*] %i/%i\n',f,c.NumTestSets);

	%%%%%%%%%%%%%%%%%%
	% ACTIVE LEARNING
	%%%%%%%%%%%%%%%%%%

	%%%
	%random initial (10%)
	%%
	avalableIdx  = find(c.training(f)==1)';			%all in the training set are avalable to be introduced to AL traning
	%removeIdxIdx = randsample( length(avalableIdx), ...
	%		        numInstStart);			%pick 10% to remove (move from avalable to AL train)
	%trainIdx     = avalableIdx(removeIdxIdx); 		%add the indexes to trainIdx
	%avalableIdx(removeIdxIdx)=[]; 				%remove avalability
	trainIdx    = [];


	%%
	% for eadch active learner
	% progress trainIdx and avalableIdx differently
	% get classification stats on the pick 
	%%
	for alType=usemethods
	    alType=alType{1};
	    fprintf('\t[%s]\n',alType);
	    
	    %initial step has no selection
	    i=0;
	    AL.(alType).ontrain.avalableIdx=avalableIdx;

	    %%%%%%
	    % start AL iterations
	    %%%%%%

	    %%%% contine until there are no points left
	    while(length(AL.(alType).ontrain.avalableIdx)>BatchSize-1)

		%%%%Instance selection for next training round
		i=i+1;
		if(i==1)
		    AL.(alType).ontrain.trainIdx=[];
		    prevLast=length(AL.(alType).ontrain.trainIdx);
		    %always use random for first
		    %AL.(alType).ontrain = updateTrainIdxRND(AL.(alType).ontrain);
		    %use what the structure says to
		    AL.(alType).ontrain = AL.(alType).initup(AL.(alType).ontrain);
		else
		    prevLast=length(AL.(alType).ontrain.trainIdx);
		    %pick what to look at next
		    AL.(alType).ontrain = AL.(alType).updater(AL.(alType).ontrain);
		end

		%what was added
		AL.(alType).added{f}{i} = AL.(alType).ontrain.trainIdx(prevLast+1:end) ;


		%how much of the last pick is postive
		AL.(alType).numPos{f}(i) = howManyPos(i,AL.(alType).ontrain);
		fprintf('[*] %i\t%i\t%i P\t',i, ...
		       length(AL.(alType).ontrain.trainIdx), ...
		       AL.(alType).numPos{f}(i) );

		%how much did it cost to make this selection
		AL.(alType).cost{f}(i) = howExpensive(i,AL.(alType).added{f}{i});
		fprintf('$%.2f', AL.(alType).cost{f}(i) );

		%%%%%MODEL
		%populate model (global) 
		AL.(alType).trainer(f,i,AL.(alType).ontrain.trainIdx);
		fprintf(' - trained',i);
		%classify this model on what is left out for testing
		AL.(alType).ontest  = AL.(alType).testCly(f,i, AL.(alType).ontest);
		fprintf(' - tested');
		%classify this model on what we've seen so far
		AL.(alType).ontrain = AL.(alType).trainCly(f,i,AL.(alType).ontrain);
		fprintf(' - evaluated\n');


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
	%AL.(alType).numInstStart=numInstStart;
	filename=[ ResultsFolder '/' whichiteration (alType) '.mat'];
	fprintf('{*} Saving to %s\n', filename);
	save(filename,'-struct', 'AL', alType);

	%%%% CHECK OUTPUT
	%fprintf('second training set looks like:\n');
	%reshape(AL.(alType).training{2}([1:10;end-9:end]),10,2)
	%fprintf('first TPs look like');
	%reshape(AL.(alType).ontrain.TP{1}([1:10;end-9:end]),10,2)

    end
    %save partion stuff
    %save('partion.mat','c');





end

%%%Forest Train %%%%%
function trainRF(f,i,trainIdxs)
    clear model
    globals;
    model = classRF_train( genes.feats(trainIdxs,:), genes.labels(trainIdxs,:),numTrees );
end

%%%SMO Train %%%%%
function trainSMO(f,i,trainIdxs)
    clear model; globals;
    model = svmtrain( genes.feats(trainIdxs,:), genes.labels(trainIdxs,:),'kernel_function', 'rbf');
end


%%% SMO test results %%%%%
function ontest = classifySMO_ontest(f,i,ontest)
    globals;
    %get total results
    results=zeros(1,c.TestSize(f));
    Idx=c.test(f);
    results = svmclassify(model,genes.feats(Idx,:));

    ontest=setResults(f,i,results,ontest,Idx);
    %fprintf('[** test]\n');
end

%%% SMO train results %%%%%
function ontrain = classifySMO_ontrain(f,i,ontrain)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    Idx=[ontrain.trainIdx ontrain.avalableIdx];
    results = svmclassify(model,genes.feats(Idx,:));

    %set TP,FP,TN,TP for each batch
    ontrain=setResults(f,i,results,ontrain,Idx);
    %fprintf('[** train] %i - %i\n', i, length(ontrain.trainIdx));
end

%%% Forset test results %%%%%
function ontest=classifyRF_ontest(f,i,ontest)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    Idx=c.test(f);
    results = classRF_predict(genes.feats(Idx,:),model);

    ontest=setResults(f,i,results,ontest,Idx);
    %fprintf('[** test] %i\n', i);
end

%%% Forset train results %%%%%
function ontrain = classifyRF_ontrain(f,i,ontrain)
    globals;

    %get total results
    results=zeros(1,c.TestSize(f));
    Idx=[ontrain.trainIdx ontrain.avalableIdx];
    results = classRF_predict(genes.feats(Idx,:),model);

    ontrain=setResults(f,i,results,ontrain,Idx);
    %fprintf('[** train] %i - %i\n', i, length(ontrain.trainIdx));
end

%%% Generic Results from Interation i on partition f given results and an index for labels
function rStruct=setResults(f,i,results,rStruct,Idx)
    %For accuracy/recall per iteration
    globals;
    posIdx = find(genes.labels(Idx) == 1);
    negIdx = find(genes.labels(Idx) == 0);

    rStruct.TP{f}(i) = length(  find( results(posIdx) ==1 )  );
    rStruct.TN{f}(i) = length(  find( results(negIdx) ==0 )  );
    rStruct.FP{f}(i) = length(  find( results(negIdx) ==1 )  );
    rStruct.FN{f}(i) = length(  find( results(posIdx) ==0 )  );
end

%%%%%%%%%%%%%%%%
% how many positives were just selected in the most recent batch?
% because trainIdx is global for each fold, this cannot be done afterward
% should changed trainIdx to be trianIdx{f} so it is more general
%%%%%%%%%%%%%%%%
function numPos = howManyPos(i,ontrain)
   global fold genes;
   %batch size bewteen this and last (should be BatchSize, unless KMeans or end of the line)
   if(i==1)
       b=ontrain.batch{fold}(1);
   else
       b=ontrain.batch{fold}(i-1:i);
       b=b(2)-b(1);
   end

   %note trainIdx isn't a cell, it is reset each fold
   %this won't work if done for f not on current fold
   mostRecentIdx=ontrain.trainIdx(end-b+1:end);

   %how many of what we just found are positive?
   posIdx = find(genes.labels(mostRecentIdx) == 1);

   numPos=length(posIdx);

   %quick check, does b-pos = neg
   %Neg=b-numPos;
   %NegL=length(find(genes.labels(mostRecentIdx)==0));
   %if(Neg~=NegL); fprintf('ERR - expected negatives (%i) != labeled (%i)\n',Neg, NegL); end;

end

%%%%%%%%%%
% What is the cost of the selection
%%%%%%%%%
function cost = howExpensive(i,mostRecentIdx)
   global fold genes;

   %%% gc and length determin cost
   len = genes.feats(mostRecentIdx,11);
   gc  = genes.feats(mostRecentIdx,12);

   cost=sum(costOne(len,gc));
end


%%%%%%% 
%% Update instances to look at
%%%%%%

%%% Random New Instances %%%
function ontrain = updateTrainIdxRND(ontrain)
    globals;
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 

    removeIdxIdx = randsample( length(ontrain.avalableIdx), batch );

    ontrain=updateTrainIdx(ontrain,removeIdxIdx);
end

%%% K-means density based %%%
function ontrain = updateTrainIdxKMeans(ontrain)
    globals;
    %%Clustring
    %% haven't built the cluster yet!
    %% Only finding centroids once, should they be found on only new data
    %% seems like this wouldn't add anything
    %%if(~isfield(ontrain,'clusterSorted'))
    %%end

    %%change my mind
    %easier to cluster every time, on the same data, to get different indexes
    %%clustering isn't that slow

    AllIdx=[ontrain.avalableIdx, ontrain.trainIdx];
    %D should be in the same length as c.train(f)
    %so we cluster on all avalable training data

    %cluster
    [clusters,~,~,distances]=kmeans(genes.feats(AllIdx,:),Clusters);
    %remove already trained indexes
    N=length(ontrain.avalableIdx);
    clusters=clusters(1:N); distances=distances(1:N,:);

    %for each cluster,sort and find
    removeIdxIdx=[]; %not sure what the size will be at the end
    for k=1:Clusters
	%find this cluster
	KclustAvIdx=find(clusters==k);

	%sort distances of all in this cluster
	[~,sortIdx]=sort(distances(KclustAvIdx));

	%grab the index
	avlIdxIdx = KclustAvIdx(sortIdx)';

	%decide how many to take
	ni=length(KclustAvIdx);
	si=floor(BatchSize*ni/N);

	%or only take what is left
	if(si>ni)
	    si=ni; 
	    fprintf('ERR -- desired portion too big %i>%i\n', si,ni);
	end

	if(si>0)
	    removeIdxIdx=[removeIdxIdx avlIdxIdx(1:si)];
	else
	    fprintf('ERR -- nothing in %i\n', k);
	end
	

    end
    if(length(removeIdxIdx) < BatchSize -10)
	fprintf('ERR -- BatchSize off by a lot: %i (%i)\n', length(removeIdxIdx),BatchSize);
    end

    ontrain=updateTrainIdx(ontrain,removeIdxIdx);

end

%%% Most votes (RF) become new instances %%%
%%% There is a segfault bug in here sometimes, I think it is size of removeIdxIdx end-batch-1?
function ontrain = updateTrainIdxLeastVotes(ontrain)
    globals;

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote difference
    [results,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    [~,votes]       = sort(abs(votes(:,1) - votes(:,2)),'descend'); %the bigger the difference, the better

    %set batch
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 


    %update usagage indexes
    removeIdxIdx     = votes(end-batch-1:end);

    ontrain=updateTrainIdx(ontrain,removeIdxIdx);
end

%%% Pick Most + Votes (RF) for new instances %%%
function ontrain = updateTrainIdxLeastPositiveVotes(ontrain)
    globals;

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote difference
    [results,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    [~,votes]       = sort(votes(:,1) - votes(:,2),'descend'); %positive votes are true class

    %set batch
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    removeIdxIdx     = votes(end-batch-1:end);

    ontrain=updateTrainIdx(ontrain,removeIdxIdx);
end

function ontrain = updateTrainIdxRFEntropy(ontrain)
    globals;

    %entropy function
    e=@(p,n) -1.*(p/numTrees.*log2(p/numTrees)+ n.*log2(n/numTrees));

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote entropy
    [~,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    votes=votes./numTrees;

    entropy=e(votes(:,1),votes(:,2));
    nanIDX=find(isnan(entropy));
    entropy(nanIDX)=zeros(size(nanIDX));

    [r,entropy]       = sort(entropy,'ascend');  %put lowest on top, pick from bottom
    %fprintf('\t\ttaking entropy %.2f,..,[..,%.2f, %.2f]\n', r([1 end-1:end]));

    %set batch size
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    removeIdxIdx = entropy(end-batch-1:end);

    ontrain      = updateTrainIdx(ontrain,removeIdxIdx);
end
function ontrain = updateTrainIdxRFConfusion(ontrain)
    globals;

    %hwo many models do we have saved?
    numMods=0;
    if(isfield(ontrain,'models')); numMods = length(ontrain.models); end;

    %keep at most 3 models around
    %FILO
    switch numMods
	case {2, 3}
	    ontrain.models={model ontrain.models{1:2}};
	    %now have 3 models, use this method!

	case 1
	    ontrain.models={model ontrain.models{1}};

	    %don't have enough models to be useful yet, use Kmeans
	    ontrain = updateTrainIdxKMeans(ontrain);
	    return

	case 0
	    ontrain.models={model};

	    %don't have enough models to be useful yet, use Kmeans
	    ontrain = updateTrainIdxKMeans(ontrain);
	    return
	otherwise
	    fprintf('ERR -- why is modlen %i\n',numMods);
    end

    %predict for each model stored
    % into votes cell 
    numMods=length(ontrain.models);
    for m=1:numMods;
	[~,votes{m}] = classRF_predict(genes.feats(ontrain.avalableIdx,:), ontrain.models{m});
	votes{m}=votes{m}./numTrees;
    end


    %%%%%%%%%%%%%
    %calcluate confusion
    %%%%%%%%%%%%%

    %average probability
    %(cat to cell arrays to 3rd demension, sum across dimension)/cell length
    avgP=sum(cat(3,votes{:}),3)./length(votes);

    %entropy function
    e=@(p,n,ap,an) (p.*log2(p./ap)+ n.*log2(n./an));

    %zero
    confusion=zeros(length(votes{1}),1);
    %for each entropy
    for m=1:numMods
	%NaN is a problem
	entropy=e(votes{m}(:,1),votes{m}(:,2),avgP(:,1),avgP(:,2));
	nanIDX=find(isnan(entropy));
	entropy(nanIDX)=zeros(size(nanIDX));

	confusion=confusion+entropy;

    end


    [r,confusion] = sort(confusion,'ascend'); %smallest on top, pick from bottom
    %fprintf('\t\ttaking confusion %.2f,..,[..,%.2f, %.2f]\n', r([1 end-1:end]));

    %set batch size
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    if(batch<BatchSize);fprintf('small batch: %i; conf: %i\n',batch,length(confusion));end
    removeIdxIdx = confusion(end-batch-1:end);

    %update indexes
    ontrain      = updateTrainIdx(ontrain,removeIdxIdx);
end

%%% Pick Most - Votes (RF) for new instances %%%
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
    removeIdxIdx     = votes(end-batch-1:end);
    ontrain=updateTrainIdx(ontrain,removeIdxIdx);
end

%%%Generic update function
%% take index of indexs to remove from avalable
%% remove, add to train, and put new size in 'batch'
function ontrain=updateTrainIdx(ontrain,removeIdxIdx)
    %globals; %so we can use fold
    global fold;
    %update indexes
    ontrain.trainIdx = [ontrain.trainIdx ontrain.avalableIdx(removeIdxIdx)];
    ontrain.avalableIdx(removeIdxIdx) = [];

    %add batch size
    if(~isfield(ontrain,'batch') || length(ontrain.batch)<fold)
	ontrain.batch{fold}=[length(removeIdxIdx)];
    else
	ontrain.batch{fold}=[ontrain.batch{fold} ontrain.batch{fold}(end)+length(removeIdxIdx)];
    end
end




%%LEGACY for name mess up
function updateTrainIdxMostPositiveVotes
end
function ontrain=updateTrainIdxMostVotes(ontrain)
end

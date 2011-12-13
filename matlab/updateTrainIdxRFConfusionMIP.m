%%NUMBER of models to use is hard coded at 3
%% don't know how to do < $NUM-1 in matlab switch :)

function ontrain = updateTrainIdxRFConfusionMIP(ontrain)
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
    classif=zeros(length(ontrain.avalableIdx),3);
    for m=1:numMods;
	[classif(:,m),votes{m}] = classRF_predict(genes.feats(ontrain.avalableIdx,:), ontrain.models{m});
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

    %%%MIPS part

    %classif>1is most likely positve, so say 2 of 3 trees
    classif   = sum(classif,2);
    posIdx = find(classif>1);
    w = max(confusion(posIdx)); %weight to add to all positive so they are choosen first!
    confusion(posIdx)=confusion(posIdx)+w;


    %%%%%


    [r,confusion] = sort(confusion,'ascend'); %smallest on top, pick from bottom
    %fprintf('\t\ttaking confusion %.2f,..,[..,%.2f, %.2f]\n', r([1 end-1:end]));

    %set batch size
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    if(batch<BatchSize);fprintf('small batch: %i; conf: %i\n',batch,length(confusion));end
    removeIdxIdx = confusion(end-batch:end);

    %update indexes
    ontrain      = updateTrainIdx(ontrain,removeIdxIdx);
end

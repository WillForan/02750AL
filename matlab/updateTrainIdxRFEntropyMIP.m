function ontrain = updateTrainIdxRFEntropy(ontrain)
    globals;

    %entropy function
    e=@(p,n) -1.*(p.*log2(p)+ n.*log2(n));

    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote entropy
    [~,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);
    votes=votes./numTrees;

    entropy=e(votes(:,1),votes(:,2));
    nanIDX=find(isnan(entropy));
    entropy(nanIDX)=zeros(size(nanIDX));

    %%%MIPS part
    %v1-v2>1 gives all predicted diff genes
    posIdx=find(votes(:,1)-votes(:,2)>0);
    w = max(entropy(posIdx)); %weight to add to all positive so they are choosen first!
    entropy(posIdx)=entropy(posIdx)+w;
    %%%%%

    %entropy is now the index of where in the initial results the value was 
    [r,entropy]       = sort(entropy,'ascend');  %put lowest on top, pick from bottom
    %fprintf('\t\ttaking entropy %.2f,..,[..,%.2f, %.2f]\n', r([1 end-1:end]));


    %set batch size
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain-1; end; 


    %update usagage indexes
    removeIdxIdx = entropy(end-batch:end);

    ontrain      = updateTrainIdx(ontrain,removeIdxIdx);
end

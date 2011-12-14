function ontrain = updateTrainIdxCheapestFirst(ontrain)
    globals;

    %set batch size
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 

    %get the cost of all
    cost = costOne(genes.feats(ontrain.avalableIdx,11) , ...
                  genes.feats(ontrain.avalableIdx,12));

    
    [~,costIdx]       = sort(cost); %first are best


    %update usagage indexes
    removeIdxIdx  = costIdx(1:batch);

    %%check
    %cost2 = costOne(genes.feats(ontrain.avalableIdx(removeIdxIdx,11)) , ...
    %              genes.feats(ontrain.avalableIdx(removeIdxIdx,12)));

    %[cost(removeIdxIdx) cost2] 
    %consistancy check
    %fprintf('spending %.2f',sum(cost(removeIdxIdx)));
    

    %check that we are picking the lowest!  
   % [cost2([1,end])' cost(removeIdxIdx([1,end]))' ] 

    ontrain = updateTrainIdx(ontrain,removeIdxIdx);
end

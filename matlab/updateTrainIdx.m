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

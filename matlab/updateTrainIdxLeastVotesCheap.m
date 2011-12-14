
%%% Most votes (RF) become new instances %%%
%%% There is a segfault bug in here sometimes, I think it is size of removeIdxIdx end-batch-1?
%%%!!!!!!
%%%!!!!!! Least difference does best, this is labeled wrong!
function ontrain = updateTrainIdxLeastVotesCheap(ontrain)
    globals;
    %set batch
    remain = length(ontrain.avalableIdx);
    batch  = BatchSize;
    if(remain < batch); batch=remain; end; 

    batchIdx=1:batch;


    %classify on only the stuff we haven't seen
    % and get Idx of avalableIdx ranked by vote difference
    [~,votes] = classRF_predict(genes.feats(ontrain.avalableIdx,:),model);

    votes=abs(votes(:,1) - votes(:,2));

    

    %move votes to an index
    [votes_sorted,votesIdx]       = sort(votes,'ascend'); %0....30
    %average of what  would  be remove
    vote_avg=mean(votes_sorted(batchIdx));
    %std dev of  ''
    vote_std=std(votes_sorted(batchIdx));


    %how far can we dip into others before we're picking things differing too much from what we want
    furthestIdx=find(votes_sorted>vote_avg+vote_std);

    if(furthestIdx)
	furthestIdx=furthestIdx(1);

	len=genes.feats(ontrain.avalableIdx(votesIdx),11);
	gc=genes.feats(ontrain.avalableIdx(votesIdx),12);

	%cost ordered by what we want to take first
	costOrderedByRank = costOne(len,gc);

	%what do I have to imporve upon to get a better cost
	[cost_max,max_idx] =max(costOrderedByRank(batchIdx));

	%search for costs outside of what we have but within our scope
	betterScope   = costOrderedByRank(batch:furthestIdx);
	betterCostIdx = find(betterScope<cost_max);

	while( betterCostIdx )
	    [mincost,mincost_idx] = min(betterScope);

	    offset=batch-1; %how much are betterScope indexes off by

	    %replace the batchIdx with the new better idx (which is offset by size of batch)
	    batchIdx(max_idx)=mincost_idx+offset;

	    %can't remove, that would mess up indexes. instead set as infinity
	    betterScope(mincost_idx)=Inf;

	    %try the same thing again
	    [cost_max,max_idx] =max(costOrderedByRank(batchIdx));
	    betterCostIdx = find(betterScope<cost_max);
	end


    
    

    end
    %update usagage indexes
    removeIdxIdx     = votesIdx(batchIdx);

    ontrain=updateTrainIdx(ontrain,removeIdxIdx);
end

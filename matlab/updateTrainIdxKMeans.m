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

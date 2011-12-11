%input
global genes;

%%% parameteres
global numTrees;
global BatchSize;
global validationFolds;
global Clusters;

%%% partions
global fold; %meh, sometimes it's passed as f; other times used as global fold


%results
global AL; %methods and results
global c; %partitions

%model
global model;
%global results votes;

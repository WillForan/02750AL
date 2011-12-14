%%%%%%%%%%%%%%%%%%%
%define methods
%%%%%%%%%%%%%%%%%%%
%Least Votes
AL.RFLeastVotes.ontest   = struct; AL.RFLeastVotes.ontrain = struct;
AL.RFLeastVotes.initup   = @updateTrainIdxRND;
AL.RFLeastVotes.updater  = @updateTrainIdxLeastVotes;
AL.RFLeastVotes.trainer  = @trainRF;
AL.RFLeastVotes.testCly  = @classifyRF_ontest;
AL.RFLeastVotes.trainCly = @classifyRF_ontrain;
AL.RFLeastVotes.imgfile  = 'RF-MV';
AL.RFLeastVotes.title    = ['Random Forest (' num2str(numTrees)      ...
			     ')- Least Votes selection' ];

AL.RFLeastVotesCheap             = AL.RFLeastVotes;
AL.RFLeastVotesCheap.updater     = @updateTrainIdxLeastVotesCheap;
AL.RFLeastVotesCheap.imgfile     = 'RF-MVC';
AL.RFLeastVotesCheap.title       = ['Random Forest (' num2str(numTrees)      ...
       		                       ')- Least Chepeast Votes' ];

AL.RFLeastPositiveVotes             = AL.RFLeastVotes;
AL.RFLeastPositiveVotes.updater     = @updateTrainIdxLeastPositiveVotes;
AL.RFLeastPositiveVotes.imgfile     = 'RF-MPV';
AL.RFLeastPositiveVotes.title       = ['Random Forest (' num2str(numTrees)      ...
       		                       ')- Least Positive Votes' ];

AL.RFRandom             = AL.RFLeastVotes;
AL.RFRandom.updater     = @updateTrainIdxRND;
AL.RFRandom.imgfile     = 'RF-RND';
AL.RFRandom.title       = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Random selection' ];

AL.RFKMeans            = AL.RFLeastVotes;
AL.RFKMeans.updater    = @updateTrainIdxKMeans;
AL.RFKMeans.imgfile    = 'RF-KM';
AL.RFKMeans.title      = ['Random Forest (' num2str(numTrees)      ...
       		     ')- KMeans selection' ];

AL.RFEntropy           = AL.RFLeastVotes;
AL.RFEntropy.updater   = @updateTrainIdxRFEntropy;
AL.RFEntropy.imgfile   = 'RF-E';
AL.RFEntropy.title     = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Entropy selection' ];

AL.RFEntropyMIP         = AL.RFLeastVotes;
%AL.RFEntropyMIP.initup  = @updateTrainIdxKMeans;
AL.RFEntropyMIP.updater = @updateTrainIdxRFEntropyMIP;
AL.RFEntropyMIP.imgfile = 'RF-EntropyMIP';
AL.RFEntropyMIP.title   = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Entropy MIP selection' ];

AL.RFConfusion         = AL.RFLeastVotes;
AL.RFConfusion.initup  = @updateTrainIdxKMeans;
AL.RFConfusion.updater = @updateTrainIdxRFConfusion;
%AL.RFConfusion.default = @updateTrainIdxKmeans; %does not work
AL.RFConfusion.imgfile = 'RF-Confusion';
AL.RFConfusion.title   = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Confusion selection' ];

AL.RFConfusionMIP         = AL.RFLeastVotes;
AL.RFConfusionMIP.initup  = @updateTrainIdxKMeans;
AL.RFConfusionMIP.updater = @updateTrainIdxRFConfusionMIP;
AL.RFConfusionMIP.imgfile = 'RF-ConfusionMIP';
AL.RFConfusionMIP.title   = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Confusion MIP selection' ];


AL.RFLV             = AL.RFLeastVotes;
AL.RFLV.updater     = @updateTrainIdxLeastVotes;
AL.RFLV.imgfile     = 'RF-LV';
AL.RFLV.title       = ['Random Forest (' num2str(numTrees)      ...
			     ')- Least Votes selection' ];

AL.SMORND.ontest   = struct; AL.RFLeastVotes.ontrain = struct;
AL.SMORND.initup   = @updateTrainIdxRND;
AL.SMORND.updater  = @updateTrainIdxRND;
AL.SMORND.trainer  = @trainSMO;
AL.SMORND.testCly  = @classifySMO_ontest;
AL.SMORND.trainCly = @classifySMO_ontrain;
AL.SMORND.imgfile  = 'SMO-RND';
AL.SMORND.title    = ['SMO - Random ' ];


AL.RFCheap          = AL.RFLeastVotes;
AL.RFCheap.initup   = @updateTrainIdxCheapestFirst;
AL.RFCheap.updater  = @updateTrainIdxCheapestFirst;
AL.RFCheap.imgfile  = 'RF-Cheapest';
AL.RFCheap.title    = ['Random Forest (' num2str(numTrees)      ...
       		           ')- Cheapest First' ];


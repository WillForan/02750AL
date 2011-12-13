%%%%%%%%%%%%%%%%%%%
%define methods
%%%%%%%%%%%%%%%%%%%
%Most Votes
AL.RFMostVotes.ontest   = struct; AL.RFMostVotes.ontrain = struct;
AL.RFMostVotes.initup   = @updateTrainIdxRND;
AL.RFMostVotes.updater  = @updateTrainIdxMostVotes;
AL.RFMostVotes.trainer  = @trainRF;
AL.RFMostVotes.testCly  = @classifyRF_ontest;
AL.RFMostVotes.trainCly = @classifyRF_ontrain;
AL.RFMostVotes.imgfile  = 'RF-MV';
AL.RFMostVotes.title    = ['Random Forest (' num2str(numTrees)      ...
			     ')- Most Votes selection' ];

AL.RFMostPositiveVotes             = AL.RFMostVotes;
AL.RFMostPositiveVotes.updater     = @updateTrainIdxMostPositiveVotes;
AL.RFMostPositiveVotes.imgfile     = 'RF-MPV';
AL.RFMostPositiveVotes.title       = ['Random Forest (' num2str(numTrees)      ...
       		                       ')- Most Positive Votes' ];

AL.RFRandom             = AL.RFMostVotes;
AL.RFRandom.updater     = @updateTrainIdxRND;
AL.RFRandom.imgfile     = 'RF-RND';
AL.RFRandom.title       = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Random selection' ];

AL.RFKMeans            = AL.RFMostVotes;
AL.RFKMeans.updater    = @updateTrainIdxKMeans;
AL.RFKMeans.imgfile    = 'RF-KM';
AL.RFKMeans.title      = ['Random Forest (' num2str(numTrees)      ...
       		     ')- KMeans selection' ];

AL.RFEntropy           = AL.RFMostVotes;
AL.RFEntropy.updater   = @updateTrainIdxRFEntropy;
AL.RFEntropy.imgfile   = 'RF-E';
AL.RFEntropy.title     = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Entropy selection' ];

AL.RFEntropyMIP         = AL.RFMostVotes;
AL.RFEntropyMIP.updater = @updateTrainIdxRFEntropyMIP;
AL.RFEntropyMIP.imgfile = 'RF-EntropyMIP';
AL.RFEntropyMIP.title   = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Entropy MIP selection' ];

AL.RFConfusion         = AL.RFMostVotes;
AL.RFConfusion.initup  = @updateTrainIdxKMeans;
AL.RFConfusion.updater = @updateTrainIdxRFConfusion;
%AL.RFConfusion.default = @updateTrainIdxKmeans; %does not work
AL.RFConfusion.imgfile = 'RF-Confusion';
AL.RFConfusion.title   = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Confusion selection' ];

AL.RFConfusionMIP         = AL.RFMostVotes;
AL.RFConfusionMIP.initup  = @updateTrainIdxKMeans;
AL.RFConfusionMIP.updater = @updateTrainIdxRFConfusionMIP;
AL.RFConfusionMIP.imgfile = 'RF-ConfusionMIP';
AL.RFConfusionMIP.title   = ['Random Forest (' num2str(numTrees)      ...
       		     ')- Confusion MIP selection' ];


AL.RFLV             = AL.RFMostVotes;
AL.RFLV.updater     = @updateTrainIdxLeastVotes;
AL.RFLV.imgfile     = 'RF-LV';
AL.RFLV.title       = ['Random Forest (' num2str(numTrees)      ...
			     ')- Least Votes selection' ];

AL.SMORND.ontest   = struct; AL.RFMostVotes.ontrain = struct;
AL.SMORND.initup   = @updateTrainIdxRND;
AL.SMORND.updater  = @updateTrainIdxRND;
AL.SMORND.trainer  = @trainSMO;
AL.SMORND.testCly  = @classifySMO_ontest;
AL.SMORND.trainCly = @classifySMO_ontrain;
AL.SMORND.imgfile  = 'SMO-RND';
AL.SMORND.title    = ['SMO - Random ' ];


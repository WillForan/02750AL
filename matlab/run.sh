#!/usr/bin/env bash

source runAL.sh

#run 5 times for each (which will run 5 fold cross validation)
#needs a folder to put mat files ... set as alResults/ now
 for d in  'PF' 'PPI'; do
  for m in 'RFRandom' 'RFConfusion'  'RFEntropyMIP' 'RFConfusionMIP' 'RFKMeans'  'RFMostVotes' 'RFMostPositiveVotes' 'RFEntropy' 'SMORND' ; do
     echo "===========starting $m $d"
     runAL $m $d 5;
     echo "===========finished $m $d"
   done
 done; 

#plot the two datasets
#make sure to make the directories they want (GeneResults/ and PPIResults/)
 matlab -nojvm -nodisplay -r "plotAL('PF'); plotAL('PPI');quit"


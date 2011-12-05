source wekaShortcuts.sh

time wekacli nb 	weka.classifiers.bayes.NaiveBayes
time wekacli tree 	weka.classifiers.trees.J48 -C 0.25 -M 2 -t
time wekacli rndforest 	weka.classifiers.trees.RandomForest -I 10 -K 0 -S 1
time wekacli kmean 	weka.clusterers.SimpleKMeans -N 2 -A "weka.core.ManhattanDistance -R first-last" -I 500 -S 10
time wekacli smo 	weka.classifiers.functions.SMO -C 1.0 -L 0.001 -P 1.0E-12 -N 0 -V -1 -W 1 -K "weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E 1.0"

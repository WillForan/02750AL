WEKAJAR=/usr/share/java/weka/weka.jar
ARFF=plasDB/labeled.arff

function wekacli {
    if [ -z "$(echo $2|grep weka)" ];then
	echo wekacli expects form
	echo wekacli treeoutputname weka.classifiers.trees.J48 -C 0.25 -M 2 
	return 
    fi
    
    outputName=wekaout/${1}-$(date +%F)
    shift

    echo $@ > $outputName
    java -Djava.awt.headless=true -classpath $WEKAJAR  $@ -t $ARFF | tee -a $outputName
}

function wekaInfo {
    while (($#)); do
	echo $1
	perl  -slane 'if(m/\[(.*?)\]/){$node{$1}++}  $c=$_ if(m/^Correct/); $i=$_ if(/^Inc/); END{print $c,"\n",$i,"\n",join(",",scalar(keys %node),sort keys %node)}' $1
	shift
    done 
}

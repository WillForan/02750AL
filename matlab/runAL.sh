#!/usr/bin/env bash

function runALnum {

    METHOD="$1"
    DATA="$2"
    RUN="$3"

    while true; do
	echo -e "n\n" | matlab -nodisplay -nosplash -nojvm -r "runAL({'$METHOD'},'$DATA',$RUN);exit;"


	if [ "$?" == "0" ]; then 
	 break; #We didn't crash!!!! leave this forever loop
	else
	 sleep 1 #maybe give the computer some time to think about what it's done
	fi
    done
}

function runAL {
    METHOD=$1
    DATA="$2"
    RUNS="$3"
    for run in $(seq 1 $RUNS); do
	echo "======= runALnum $METHOD $DATA $run ============"
	runALnum $METHOD $DATA $run
    done
}

function runALCLI {
    METHOD=$1
    DATA="$2"
    RUNS="$3"
    while true; do
	echo -e "n\n" | matlab -nodisplay -nojvm -r "for i=1:$RUNS; runAL({'$METHOD'},'PF',i); end;quit"
	if [ "$?" == "0" ]; then 
	 break; #We didn't crash!!!! leave this forever-loop
	else
	 sleep 1 #maybe give the computer some time to think about what it's done
	fi
    done
}



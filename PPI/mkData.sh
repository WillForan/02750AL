#!/usr/bin/env bash

wget \
http://www.cs.cmu.edu/afs/cs.cmu.edu/project/structure-9/PPI/proteins05-data-share/phyInteract/dipsPosPair.feature \
http://www.cs.cmu.edu/afs/cs.cmu.edu/project/structure-9/PPI/proteins05-data-share/phyInteract/dipsRandpairSub23w.feature

cat dipsPosPair.feature <(shuf dipsRandpairSub23w.feature|head -n 6685)|shuf > data.csv

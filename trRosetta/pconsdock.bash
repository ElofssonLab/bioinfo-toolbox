#!/bin/bash 


dockQ="/home/arnee/git/DockQ/DockQ.py"


for i in $@
do
    name=`basename $i .pdb`
    dir=`dirname $i `
    for j in $@
    do
	python3 ~/git/bioinfo-toolbox/DockQ/DockQ.py -short $i $j
    done | gawk '{if ($2<1) {i+=1;s+=$2}}END{print i,s,s/i}' > ${dir}/${name}.pconsdock
done

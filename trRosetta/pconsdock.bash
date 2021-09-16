#!/bin/bash 


dockQ="/home/arnee/git/DockQ/DockQ.py"


for i in $@
do
    name=`basename $i .pdb`
    dir=`dirname $i `
    for j in $@
    do
	python3 ~/git/DockQ/DockQ-mod.py -short $i $j
    done | gawk '{if ($2<=1) {i+=1;s+=$2}}END{ if (i>0){print i,s,s/i}else{print 0,0,0}}' > ${dir}/${name}.pcd
done

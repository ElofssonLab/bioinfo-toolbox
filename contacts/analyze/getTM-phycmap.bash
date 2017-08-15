#!/bin/bash 


TMscore=/pfs/nobackup/home/a/arnee/Software/TMscore/TMscore

for k in 20 30 40
do
    for l in 20 30 40
    do
	
	for i in  PF*/native.aligned_2.pdb
	do
	    j=`dirname $i`
	    echo -n $j "   "
	    $TMscore $i $j/$j*${k}_${l}*confold_mem/stage1/*model1.pdb | grep TM-score | grep d0
	    echo " "
	done | grep TM-score | sed s/TM-score\ =//g | sed s/\(d0=......//g > phycmap_${k}_${l}.txt
    done
done


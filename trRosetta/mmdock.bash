#!/bin/bash 

MM="/home/arnee/git/bioinfo-toolbox/trRosetta/MMalign"

for i in $@
do
    name=`basename $i .pdb`
    dir=`dirname $i `
    for j in $@
    do
	$MM $i $j | grep TM-score | head -1
    done | gawk '{if ($2<1) {i+=1;s+=$2}}END{print i,s,s/i}' > ${dir}/${name}.MMdock
done

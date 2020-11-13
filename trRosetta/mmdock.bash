#!/bin/bash 
MM="/home/arnee/git/bioinfo-toolbox/trRosetta/MMalign"

for i in $@
do
    name=`basename $i .pdb`
    dir=`dirname $i `
    for j in $@
    do
	$MM $i $j | grep TM-score | head -1
    done  | gawk 'BEGIN{s=0;i=0.00000000001}{if ($2<1) {s+=$2}{i+=1}}END{print i,s,s/i}' > ${dir}/${name}.MMdock
done

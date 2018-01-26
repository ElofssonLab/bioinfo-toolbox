#!/bin/bash 
for k in homo hetero
do
    for j in gdca 0.02.plm20
    do
	for i in HH1 HH0.0001 HH1.e-10
	do
	    echo -n $i $j $k " "
	    grep Zscore out/*$i.$j.$k.out  | gawk '{ sum+=$6;num++}; END{print num, sum/num}'
	done
    done
done

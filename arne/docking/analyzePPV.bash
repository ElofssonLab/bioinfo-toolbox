#!/bin/bash -x
for k in homo hetero
do
    for j in gdca 0.02.plm20
    do
	for i in HH1 HH0.0001 HH1.e-10
	do
	    echo $i $j $k
	    grep PPV out/*$i.$j.$k.out  | gawk '{if ($6>0.2) print $0}' |wc
	    ls out/*$i.$j.$k.out | wc
	done
    done
done

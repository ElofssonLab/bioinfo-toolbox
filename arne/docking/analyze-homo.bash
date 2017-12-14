#!/bin/bash -x

bin=/pfs/nobackup/home/a/arnee/contactpreds/docking/TM/bin
for ali in HH1 HH0.0001 HH1.e-10
do
    for dca in gdca 0.02.plm20
    do
	if [ -s seq_chains/$1-$2.${ali}.${dca} ]
	then
	    if [ ! -s out/$1-$2.${ali}.${dca}.out ]
	    then
		$bin/ppv-homodimer.py  seq_chains/$1.seq seq_chains/$2.seq tm-chains/$1.pdb.pdb-Bvalue.pdb tm-chains/$2.pdb.pdb-Bvalue.pdb seq_chains/$1.seq.${ali}.${dca} --print_dist --outfile out/$1-$2.${ali}.${dca}.dist > out/$1-$2.${ali}.${dca}.out
	    fi
	fi
    done
done
       

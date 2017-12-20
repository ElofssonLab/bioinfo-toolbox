#!/bin/bash -x

bin=/pfs/nobackup/home/a/arnee/contactpreds/docking/TM/bin
for ali in HH1 HH0.0001 HH1.e-10
do
    for dca in gdca 0.02.plm20
    do
	if [ -s seq_chains/$1-$2.${ali}.${dca} ]
	then
	    if [ ! -s out/$1-$2.${ali}.${dca}.hetero.out ]
	    then
		if [ -s seq_chains/$1.seq ]
		then
		    $bin/ppv-heterodimer.py  seq_chains/$1.seq seq_chains/$2.seq tm-chains/$1.pdb.pdb-Bvalue.pdb tm-chains/$2.pdb.pdb-Bvalue.pdb seq_chains/$1-$2.${ali}.${dca} --print_dist --outfile out/$1-$2.${ali}.${dca}.hetero.dist > out/$1-$2.${ali}.${dca}.hetero.out
		else
		    $bin/ppv-heterodimer.py  seq/$1.fa seq/$2.fa pdb_chains/$1.pdb.pdb-Bvalue.pdb tm-chains/$2.pdb.pdb-Bvalue.pdb seq/$1-$2.${ali}.${dca} --print_dist --outfile out/$1-$2.${ali}.${dca}.hetero.dist > out/$1-$2.${ali}.${dca}.hetero.out
		fi
	    fi
	fi
    done
done
       

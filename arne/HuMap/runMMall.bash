#!/bin/bash -x


name=$1
code=$2

if [ ! -f pdb/$1/$1.MMall2 ]
then
    touch -f pdb/$1/$1.MMall2 
    for i in cif/$2_*.pdb # cif/$2.cif
    do
	j=`basename $i .pdb `
	if [ ! -f pdb/$1/$1_$j.MMall ]
	then
	    touch pdb/$1/$1_$j.MMall
	    bin/runMMscore.bash pdb/$1/$1.pdb $i > pdb/$1/$1_$j.MMall
	fi
    done

    # Find Max
    grep -H TM-score= pdb/$1/${name}*.MM* | sort -nk 2 | tail -1 |  sed "s/.*\///g" | sed "s/.MM.*:/ /g" | sed "s/\_/ /g" > pdb/$1/$1.MMall2
fi

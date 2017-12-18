#!/bin/bash -x
#SBATCH -A SNIC2017-11-7
#SBATCH --output=hh.%A_%a.out
#SBATCH --error=hh.%A_%a.out
#SBATCH --array=1-520
#SBATCH -c 6
#SBATCH -t 06:00:00
#SBATCH -A SNIC2017-11-7

offset=$2
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos list.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

#id="pdbseq/1a3aB.fa"

CPU=6
directory=`dirname $id`
dir=`pwd`/$directory
foo=`basename $id .fa`
idname=`basename $foo .seq`


HH=/pfs/nobackup/home/a/arnee/git/PconsC3/extra/arne/MSA/runhhblits.py
a3m=/pfs/nobackup/home/a/arnee/git/PconsC3/extra/arne/MSA/a3mToTrimmed.py

for E in 1  0.0001 1.e-10
do
    file=$id.HH$E.a3m
    if [ ! -s ${file}  ]
    then
	${HH} -c $CPU -e ${E} -name HH${E} $id
	${a3m} -o -name $idname $id.HH${E}.a3m > $id.HH${E}.trimmed
    fi

done

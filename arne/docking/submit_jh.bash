#!/bin/bash -x
#SBATCH -A SNIC2017-11-7
#SBATCH --output=jh.%A_%a.out
#SBATCH --error=jh.%A_%a.out
#SBATCH --array=1-681
#SBATCH -c 6
#SBATCH -t 06:00:00
#SBATCH -A SNIC2017-11-7
#SBATCH -p largemem

offset=$2
list=$1

if [ ! $offset ]
then
    offset=0
fi

pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos list.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

#id="seq_chains/1a0s_P.seq"

CPU=6
directory=`dirname $id`
dir=`pwd`/$directory
foo=`basename $id .fa`
idname=`basename $foo .seq`


JH=/pfs/nobackup/home/a/arnee/git/PconsC3/extra/arne/MSA/runjackhmmer.py
a3m=/pfs/nobackup/home/a/arnee/git/PconsC3/extra/arne/MSA/a3mToTrimmed.py

for E in 1  0.0001 1.e-10
do
    file=$id.JH$E.a3m
    if [ ! -s ${file}  ]
    then
	${JH} -c $CPU -e ${E} -name JH${E} $id -db /pfs/nobackup/home/a/arnee/data/uniref100.fasta
	${a3m} -o -name $idname $id.JH${E}.a3m > $id.JH${E}.trimmed
    fi

done

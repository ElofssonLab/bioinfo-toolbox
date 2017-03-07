#!/bin/bash -x
#SBATCH --output=confold.slurm.out
#SBATCH --error=confold.slurm.error
#SBATCH -n 1
#SBATCH -t 1-00:00:00

ulimit -s unlimited

#id=$1
#echo $id

#confold=/proj/bioinfo/software/confold/CONFOLD/confold.pl
confold=$HOME/git/bioinfo-toolbox/arne/CONFOLD/confold_mem.pl
l=2.5

#dir=data/29.0/$id

#seq=$1dir/${id}.fa
#rr=$dir/${id}.hhE0.pconsc3.l3.rr
#ss=${seq}.confold.ss2

id=`basename $1 .seq`
seq=$1
rr=$2
ss=$3
top=$4
outdir=${rr}_confold_mem

if [ ! -f $outdir/stage2/${id}.fa_model1.pdb ]; then
    rm -rf $outdir
    mkdir $outdir
    $confold -top $top -seq $seq -rr $rr -ss $ss -o $outdir -selectrr ${l}L -stage2 1 -mcount 20 &> $outdir/${id}.log
fi
ls -l $outdir/stage2/${id}.fa_model1.pdb

wait

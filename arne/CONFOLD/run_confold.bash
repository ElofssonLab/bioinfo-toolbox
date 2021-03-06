#!/bin/bash -x
#SBATCH --output=confold.slurm.out
#SBATCH --error=confold.slurm.error
#SBATCH -n 1
#SBATCH -t 1-00:00:00

ulimit -s unlimited

#id=$1
#echo $id

if [ -s /proj/bioinfo/software/confold/CONFOLD/confold_mem.pl ]
then
    confold=/proj/bioinfo/software/confold/CONFOLD/confold.pl
elif  [ -s /pfs/nobackup/home/m/mircomic/confold/CONFOLD/confold.pl ]
then
    confold=/pfs/nobackup/home/m/mircomic/confold/CONFOLD/confold.pl
elif  [ -s $HOME/git/bioinfo-toolbox/CONFOLD/confold.pl ]
then
    confold=$HOME/git/bioinfo-toolbox/CONFOLD/confold.pl 
fi

#l=2.5

#dir=data/29.0/$id

#seq=$1dir/${id}.fa
#rr=$dir/${id}.hhE0.pconsc3.l3.rr
#ss=${seq}.confold.ss2

id=`basename $1 .seq | sed s/.fa$//g`
seq=$1
rr=$2
ss=$3
cutoff=$4
outdir=${rr}_confold
out=`echo ${rr} | sed s/^.*.fa.//g | sed s/.pconsc3//g`
outdir=${out}_${cutoff}_cm

if [ ! -f $outdir/stage1/${id}.fa_model1.pdb ]; then
    rm -rf $outdir
    mkdir $outdir
    $confold -seq $seq -rr $rr -ss $ss -o $outdir -selectrr ${cutoff}L -stage2 0 -mcount 50 &> $outdir/${id}.log
fi
ls -l $outdir/stage1/${id}.fa_model1.pdb

wait

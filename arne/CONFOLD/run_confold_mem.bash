#!/bin/bash -x
#SBATCH --output=confold.slurm.out
#SBATCH --error=confold.slurm.error
#SBATCH -n 1
#SBATCH -t 1-00:00:00

ulimit -s unlimited

#id=$1
#echo $id

#confold=/proj/bioinfo/software/confold/CONFOLD/confold.pl
<<<<<<< HEAD
confold=$HOME/git/bioinfo-toolbox/arne/CONFOLD/confold_mem.pl
l=2.5
=======
#confold=$HOME/git/bioinfo-toolbox/arne/CONFOLD/confold_mem.pl

if [ -s /proj/bioinfo/software/confold/CONFOLD/confold_mem.pl ]
then
    confold=/proj/bioinfo/software/confold/CONFOLD/confold_mem.pl
elif  [ -s /pfs/nobackup/home/m/mircomic/confold/CONFOLD/confold_mem.pl ]
then
    confold=/pfs/nobackup/home/m/mircomic/confold/CONFOLD/confold_mem.pl
else
    confold=$HOME/git/bioinfo-toolbox/arne/CONFOLD/confold_mem.pl
fi



#l=2.5
>>>>>>> 902382769cd14246316c6099f71f69d05967cb08

#dir=data/29.0/$id

#seq=$1dir/${id}.fa
#rr=$dir/${id}.hhE0.pconsc3.l3.rr
#ss=${seq}.confold.ss2

id=`basename $1 .seq`
seq=$1
rr=$2
ss=$3
top=$4
cutoff=$5
mindist=$6
maxdist=$7

outdir=${rr}_${cutoff}_${mindist}_${maxdist}_confold_mem

if [ ! -f $outdir/stage2/${id}.fa_model1.pdb ]; then
    rm -rf $outdir
    mkdir $outdir
    $confold -top $top -seq $seq -rr $rr -ss $ss -o $outdir -selectrr $cutoff -stage2 0 -mcount 50 -mindist $mindist -maxdist $maxdist &> $outdir/${id}.log
fi
ls -l $outdir/stage2/${id}.fa_model1.pdb

wait

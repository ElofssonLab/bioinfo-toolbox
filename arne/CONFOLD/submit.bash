#!/bin/bash -x
#SBATCH -A snic2016-10-22
# We actually start 6 jobs in parallel.
# Probably more efficient than running with 5 thread.
#SBATCH -n 6
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH -J RunJackhmmer
#SBATCH --output=out/PconsC3.%J.out
#SBATCH --error=err/PconsC3.%J.err

snic=snic2016-10-22
minitime="04:00:00"
shorttime="12:00:00"
longtime="24:00:00"
mem="2GB"

l=2.5
m=20
n=40

for i in "$@"
do
    cd $i
    j=`ls $i*.fa`
    k=`basename $j .fa`

    echo $i $j $k
    if [ ! -s ${k}.rr_${l}_${n}_${m}_confold_mem/. ]
    then
	srun --mem $mem -A $snic --time=$minitime -n 1 -c 1  $HOME/git/bioinfo-toolbox/arne/CONFOLD/run_confold_mem.bash  $k.fa $k.rr $k.confold.ss $k.top $l $n $m   &> run_confold_mem_${l}_${n}_${m}.out &
	sleep 10
    fi
    cd ..
done

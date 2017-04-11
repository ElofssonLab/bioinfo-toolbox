#!/bin/bash -lx
#SBATCH --output=parse.%A_%a.out
#SBATCH --error=parse.%A_%a.out
#SBATCH --array=1-335
#SBATCH -c 1
#SBATCH -t 30:00
#SBATCH -A SNIC2016-10-22


offset=$2
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

id="PF00001.18"

dir=`pwd`/$id

scratch=$SNIC_TMP/arnee/Summary/$id/
mkdir -p $scratch
cd $scratch


sleep 2 # waiting for filesystem

for i in $dir/*cm.tar.gz $dir/conf*[04].tar.gz $dir/*_mem.tar.gz
do
    j=`basename $i .tar.gz`
    if [ ! -s $dir/${j}_summary.out ]
    then
	if [ -s $dir/${j}_proq3.tar.gz ]
	then
	    tar -zxf $dir/${j}_proq3.tar.gz
	    ln -fs $dir/$j*.out ./
	    sleep 1
	    $dir/../bin/GetAllScores.py ${j}.tar.gz >  $dir/${j}_summary.csv
	    #	mv ${j}*out $dir/
	    sleep 1
	    rm -r ${j}/
	fi
    fi
done

cd $dir/../


#!/bin/bash -lx
#SBATCH --output=pcons.%A_%a.out
#SBATCH --error=pcons.%A_%a.out
#SBATCH --array=1-335
#SBATCH -c 1
#SBATCH -t 30:00
#SBATCH -A SNIC2016-10-22

offset=$2
list=$1
if [ -z $offset ]
then
    offset=0
fi
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

#id="PF00001.18"

dir=`pwd`/$id

scratch=$SNIC_TMP/arnee/Pcons/$id/
mkdir -p $scratch
cd $scratch


sleep 2 # waiting for filesystem

if [ ! -s $dir/${id}.raw ]
then
    rm -f qa.input
    for i in $dir/*cm.tar.gz # $dir/conf*[04].tar.gz # $dir/*_mem.tar.gz
    do
	j=`basename $i .tar.gz`
	tar -zxf $i
	ls $j/stage1/${id}*fa_[0-9].pdb >> qa.input
	ls $j/stage1/${id}*fa_[0-9][0-9].pdb >> qa.input
    done
    if [ -s qa.input ]
    then
	/pfs/nobackup/home/a/arnee/git/Pcons/bin/pconsAE.LGscore.b-an01.hpc2n.umu.se -i ./qa.input -A > ${dir}/${id}.raw
	sleep 2s
	# We do eed to split this file I think
	
    fi
    rm -r ${j}*/
fi

cd $dir/../


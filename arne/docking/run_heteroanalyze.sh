#!/bin/bash -lx
#SBATCH --output=analyze.%A_%a.out
#SBATCH --error=analyze.%A_%a.out
#SBATCH --array=1-223
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -A SNIC2017-11-7

list=$1
offset=$2

pos=$(($SLURM_ARRAY_TASK_ID + $offset))

pair=`tail -n+$pos $list | head -n1`
#pair="1a0s_P 1a0s_Q"
#cd $scratch

dir=`pwd`
scratch=$SNIC_TMP/arnee/anal/$dir/
#mkdir -p $scratch
#cd $dir


$dir/bin/analyze-hetero.bash $pair

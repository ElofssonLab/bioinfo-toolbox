#!/bin/bash -lx
#SBATCH --output=pcons.%A_%a.out
#SBATCH --error=pcons.%A_%a.out
#SBATCH --array=1-1000
#SBATCH -c 1
#SBATCH -t 30:00
#SBATCH -A SNIC2016-10-22

offset=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos tmp | head -n1`
#id=`tail -n+$pos IDs_29.0_nopdb_0.75_E3.txt | head -n1`
#id=`tail -n+$pos IDs_29.0_nopdb_0.75_E3.txt | head -n1`

offset=$2
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

id="PF00001.18"

dir=`pwd`/$id

scratch=$SNIC_TMP/arnee/TM/$id/
mkdir -p $scratch
cd $scratch


sleep 2 # waiting for filesystem

#!/bin/bash 


TMscore=/pfs/nobackup/home/a/arnee/Software/TMscore/TMscore

for i in  PF*/native.aligned_2.pdb
do
    j=`dirname $i`
    echo -n $j "   "
    $TMscore $i $j/confold_2.5_m50_hhE0/stage1/${j}*model1*stage1.pdb | grep TM-score | grep d0 | sed s/TM-score\ =//g | sed s/\(d0=......//g 
    echo " "
done


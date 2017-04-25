#!/bin/bash -lx
#SBATCH --output=TM.%A_%a.out
#SBATCH --error=TM.%A_%a.out
#SBATCH --array=1-335
#SBATCH -c 1
#SBATCH -t 90:00
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
#id="PF00002.21"
#id="PF06808.9"

dir=`pwd`/$id

scratch=$SNIC_TMP/arnee/TM/$id/
mkdir -p $scratch
cd $scratch


sleep 2 # waiting for filesystem




TMscore=/pfs/nobackup/home/a/arnee/Software/TMscore/TMscore

for i in $dir/*cm.tar.gz $dir/conf*[04].tar.gz # $dir/*_mem.tar.gz
do
    j=`basename $i .tar.gz`
    l=`grep -c PF $dir/${j}_TM.out `
#    if [ ! -s $dir/${j}_TM.out ]
    if [ $l -lt 10 ]
    then
	echo $j
	tar -zxf $i
	cp $dir/native.aligned_2.pdb .
	sleep 1
	for k in $j/stage1/${id}*.pdb
	do
	    if [ -s $k ]
	    then
		echo $k
		echo -n $k "   "  >> ${j}_TM.out
		$TMscore native.aligned_2.pdb $k | grep TM-score | grep d0 | sed s/TM-score\ =//g | sed s/\(d0=......//g  >> ${j}_TM.out
	    fi
	done
	if [ -s ${j}_TM.out ]
	then
	    cp ${j}_TM.out  $dir/${j}_TM.out  
	    sleep 2
	fi
	rm -r ${j}/
    fi
done


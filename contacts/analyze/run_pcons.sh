#!/bin/bash -lx
#SBATCH --output=pcons.%A_%a.out
#SBATCH --error=pcons.%A_%a.out
#SBATCH --array=1-335
#SBATCH -c 1
#SBATCH -t 30:00
#SBATCH -A SNIC2016-10-22

offset=$2
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

#id="PF00001.18"

dir=`pwd`/$id

scratch=$SNIC_TMP/arnee/Pcons/$id/
mkdir -p $scratch
cd $scratch


sleep 2 # waiting for filesystem

for i in $dir/*cm.tar.gz $dir/conf*[04].tar.gz $dir/*_mem.tar.gz
do
    j=`basename $i .tar.gz`
    if [ ! -s $dir/${j}_pcons.out ]
    then
	tar -zxf $i
	echo "running"
	ls $j/stage1/${id}*.pdb > qa.input
	/pfs/nobackup/home/m/mircomic/Pcons/bin/pcons -i ./qa.input -A > ${dir}/${j}.raw
#	python $dir/../bin/parse_pcons.py $j.raw > ${j}_local.out
	python $dir/../bin/reformat_pcons_out.py ${dir}/${j}.raw > ${dir}/${j}_pcons.out
#	mv ${j}*out $dir/
	sleep 2
	rm -r ${j}/
    fi
done

cd $dir/../


#!/bin/bash -lx
#SBATCH --output=pcons.%A_%a.out
#SBATCH --error=pcons.%A_%a.out
#SBATCH --array=1-335
#SBATCH -c 1
#SBATCH -t 60:00
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

for i in $dir/*cm.tar.gz # $dir/conf*[04].tar.gz # $dir/*_mem.tar.gz
do
    j=`basename $i .tar.gz`
    if [ ! -s $dir/${j}.raw ]
    then
	tar -zxf $i
	echo "running"
	ls $j/stage1/${id}*fa_[0-9].pdb > qa.input
	ls $j/stage1/${id}*fa_[0-9][0-9].pdb >> qa.input
	if [ -s qa.input ]
	then
	    /pfs/nobackup/home/m/mircomic/Pcons/bin/pcons -i ./qa.input -A > ${dir}/${j}.raw
	    #	python $dir/../bin/parse_pcons.py $j.raw > ${j}_local.out
	    # python $dir/../bin/reformat_pcons_out.py ${dir}/${j}.raw > ${dir}/${j}_pcons.out
	    #	mv ${j}*out $dir/
	    sleep 2
	fi
	rm -r ${j}/
    fi
done

cd $dir/../


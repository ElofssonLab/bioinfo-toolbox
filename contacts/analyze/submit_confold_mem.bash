#!/bin/bash -x
#SBATCH -A SNIC2016-10-22
#SBATCH --output=confold.%A_%a.out
#SBATCH --error=confold.%A_%a.out
#SBATCH --array=1-2
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -A SNIC2016-10-22

snic=snic2016-10-22
minitime="04:00:00"
shorttime="12:00:00"
longtime="24:00:00"
mem="2GB"

#offset=$2
offset=0  # Maximum  umber of 
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos $list | head -n1`

#id="PF00001.18"
#id="PF01770.15"
#id=PF00230.17
#id="PF00115.17"

dir=`pwd`/$id

#scratch=$SNIC_TMP/arnee/Summary/$id/
#mkdir -p $scratch
#cd $scratch



l=2.5
m=20
n=40

cd $dir

j=`ls ${id}*.fa`
k=`basename $j .fa`
top=${k}.fa.hhE0.pconsc3.l3.rr
if [ -s $k.confold.ss ]
then
    ss=$k.confold.ss
elif [ -s $j.confold.ss ]
then
    ss=$j.confold.ss
else
    break
fi
out=`echo ${top} | sed s/^.*.fa.//g | sed s/.pconsc3//g`
outdir=${out}_${l}_cm
	

if [ ! -s  ${outdir}.tar.gz ]
   then
       /pfs/nobackup/home/a/arnee/git/bioinfo-toolbox/arne/CONFOLD/run_confold.bash  $k.fa $top ${ss}  $l   &> ${outdir}.out
       sleep 2
       if [  -s ${outdir}/stage1/${k}.fa_1.pdb  ]
       then
	   tar -zcvf ${outdir}.tar.gz ${outdir}  --remove-files		
       fi
       sleep 2
fi

#exit 0

for m in 20 30 40 50
do
    for n in 20 30 40 50
    do
	echo ${id} $j $k
	
	out=`echo ${top} | sed s/^.*.fa.//g | sed s/.pconsc3//g`
	outdir=${out}_${l}_${n}_${m}_cm
	
	
	if [ -s ${outdir}.tar.gz ]
	then
	    tar -zxvf ${outdir}.tar.gz 
	fi
	if [ ! -s ${outdir}/stage1/${k}.fa_1.pdb  ]
	then
	    if [ ! -s  $k.fa ]
	    then
		break
	    fi
	    if [ ! -s $top ]
	    then
		break
	    fi
	    if [ ! -s   ${ss} ]
	    then
		break
	    fi
	    if [ ! -s $k.top ]
	    then
		break
	    fi
	    touch  ${outdir}/stage1/${k}.fa_1.pdb
	    /pfs/nobackup/home/a/arnee/git/bioinfo-toolbox/arne/CONFOLD/run_confold_mem.bash  $k.fa $top ${ss} $k.top $l $n $m   &> ${outdir}.out
	    sleep 2
	    if [  -s ${outdir}/stage1/${k}.fa_1.pdb  ]
	    then
		tar -zcvf ${outdir}.tar.gz ${outdir}  --remove-files		
	    fi
	    sleep 2
	fi
    done
done

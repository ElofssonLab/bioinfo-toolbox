#!/bin/bash -x
#SBATCH -A SNIC2017-11-7
#SBATCH --output=confold.%A_%a.out
#SBATCH --error=confold.%A_%a.out
#SBATCH --array=1-1
#SBATCH -c 1
#SBATCH -t 12:00:00

snic=SNIC2017-11-7
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

for i in CASP_PconsC4/P4_${id}/*p4.trimmed
do
    j=`basename $i .p4.trimmed`
    k=$id
    if [ ! -s $k/$j.rr ]
    then
	bin2/addheader.bash $k/$k.seq $i > $k/$j.rr
    fi
done


l=2.5

cd $dir

#j=`ls ${id}/*.seq`
#k=`basename $j .seq`
j=${id}.seq
k=${id}

for i in *.rr # *p4
do 
    ss=${j}.JH0.001.ss
    if [ ! -s $ss ]
    then
	/home/a/arnee/git/bioinfo-toolbox/arne/psipred_to_ss.pl ${ss}2 > ${ss}
    fi
    outdir=${k}_${i}_${l}
	

    if [ ! -s  ${outdir}.tar.gz ]
    then
	/pfs/nobackup/home/a/arnee/git/bioinfo-toolbox/arne/CONFOLD/run_confold.bash  $j $i  ${ss}  $l   &> ${outdir}.out
	sleep 2
	if [  -s ${outdir}/stage1/${k}.fa_1.pdb  ]
	then
	    echo "tar -zcvf ${outdir}.tar.gz ${outdir}  --remove-files		"
	fi
	sleep 2
    fi

    #exit 0
    
done
    

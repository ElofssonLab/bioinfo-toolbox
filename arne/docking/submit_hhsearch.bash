#!/bin/bash -x
#SBATCH -A SNIC2017-11-7
#SBATCH --output=hhsearch.%A_%a.out
#SBATCH --error=hhsearc.%A_%a.out
#SBATCH --array=1-600
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -A SNIC2017-11-7

#offset=$2
offset=0  # Maximum  umber of 
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos list.txt | head -n1`
id=`tail -n+$pos $list | head -n1`


#id="./sbf36_2_tm_rep/I2EX85/I2EX85.fasta"

directory=`dirname $id`
dir=`pwd`/$directory
idname=`basename $id .fa`

#scratch=$SNIC_TMP/arnee/hhsearch/${idname}/
#mkdir -p $scratch

#cd $scratch

#cd $dir

HHSEARCH=/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/bin/hhsearch
DB=/pfs/nobackup/home/a/arnee/data/hhsearch_dbs/pdb70_mmcif_06sep17/pdb70_hhm_db

for MSA in HH1  HH0.0001 HH1.e-10
do
    a3m=$id.$MSA.a3m
    hhr=$id.$MSA.hhr
    if [ ! -s ${hhr}  ]
    then
	${HHSEARCH} -i $a3m -d $DB
    fi

done

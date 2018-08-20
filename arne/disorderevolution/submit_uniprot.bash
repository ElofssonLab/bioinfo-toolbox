#!/bin/bash -x
#SBATCH -A SNIC2018-1-26
#SBATCH --output=uniprot.%A_%a.out
#SBATCH --error=uniprot.%A_%a.out
#SBATCH --array=1-1000
#SBATCH -c 1
#SBATCH -t 24:00:00

snic=SNIC2018-1-26
time="24:00:00"
mem="2GB"

offset=$2
#offset=0  # Maximum  number of 
list=$1
pos=$(($SLURM_ARRAY_TASK_ID + $offset))
#id=`tail -n+$pos IDs_29.0_test_done_300.txt | head -n1`
id=`tail -n+$pos $list | head -n1`
#id="./data/proteomes/UP000196581_946573.fasta"

bar=/scratch/bar.arnee.$$
foo=/scratch/foo.arnee.$$

dir="/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/bin/"

outdir="/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/data/proteomes/"

j=`basename ${id} .fasta`
if [ ! -s ${outdir}/${j}.txt ]
then
    for k in `grep \> ${id} | sed "s/|/ /g" | gawk '{print $2}'`
    do
	wget -qO - https://www.uniprot.org/uniprot/${k}.txt  >>  ${outdir}/$j.txt
    done
else
    # This was too slow.
    #    for k in `grep \> ${id} | sed "s/|/ /g" | gawk '{print $2}'`
    grep  \> $id | gawk '{print $1}' | sed  "s/.*|//g"  | sort > ${foo}
    grep -E ^ID ${outdir}/$j.txt | gawk '{print $2}' | sort > ${bar}
    for k in `diff ${foo} ${bar} | grep \< | gawk '{print $2}'`
    do
	wget -qO - https://www.uniprot.org/uniprot/${k}.txt  >>  ${outdir}/$j.txt
    done
fi

#if [ ! -s ${outdir}/${j}.xml ]
#then
#    for k in `grep \> ${id} | sed "s/|/ /g" | gawk '{print $2}'`
#    do
#	wget -qO - https://www.uniprot.org/uniprot/${k}.xml  >>  ${outdir}/$j.xml
#    done
#fi

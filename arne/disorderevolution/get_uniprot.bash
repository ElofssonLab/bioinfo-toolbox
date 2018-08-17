#!/bin/bash -x


# Just takes a fasta file of a proteome as inpu

dir=/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/data/proteomes

j=`basename $1 .fasta`
if [ ! -s $j.xml ]
then
    for k in `grep \> $i | sed "s/|/ /g" | gawk '{print $2}'`
    do
	wget -qO - https://www.uniprot.org/uniprot/${k}.xml  >>  ${dir}/$j.xml
    done
done

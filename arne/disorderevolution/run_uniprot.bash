#!/bin/bash 

for id in `cat $1`
do
    bar=/tmp/bar.arnee.$$
    foo=/tmp/foo.arnee.$$

    dir="/scratch2/arne/annotate_uniprot_proteomes/bin/"

    outdir="/scratch2/arne/annotate_uniprot_proteomes/data/proteomes/"

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
	    echo $j,$k
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
done

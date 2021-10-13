#!/bin/bash -x

while I= read -r i j k l n o
do
    if [ ! -f pdb/${i}/${i}.DockQall2 ]
    then
	touch pdb/${i}/${i}.DockQall2
	m=`basename $j .pdb`
	bin/matchpdb.py pdb/${i}/${i}.pdb pdb/${i}/${m}.pdb
	cp  pdb/${i}/${m}.pdb pdb/${i}/${m}_extra.pdb
	bin/matchpdb.py pdb/${i}/${i}_reorder.pdb pdb/${i}/${m}_extra.pdb
	pdb_reres pdb/${i}/${m}_matching.pdb  > pdb/${i}/${m}_renum.pdb 
	pdb_reres pdb/${i}/${m}_extra_matching.pdb > pdb/${i}/${m}_extra_matching_renum.pdb
	pdb_reres pdb/${i}/${i}_matching.pdb >pdb/${i}/${i}_matching_renum.pdb
	pdb_reres pdb/${i}/${i}_reorder_matching.pdb >pdb/${i}/${i}_reorder_matching_renum.pdb
	#pdb_reres $2 > $B
	#
	if [ -f pdb/${i}/${i}*.MMall2 ]
	then
	    for p in pdb/${i}/${i}*.MMall2
	    do
		pdb=`gawk '{print $2}' $p `
		chain=`gawk '{print $3}' $p `
		if [ -f cif/${pdb}_${chain}_rechain.pdb ]
		then
		    cp cif/${pdb}_${chain}_rechain.pdb pdb/${i}/${m}_mmall.pdb 
		fi
	    done
	fi
	if [ -f cif/${k}_${l}${n}_rechain.pdb ]
	then
	    cp cif/${k}_${l}${n}_rechain.pdb pdb/${i}/${m}_pdborg.pdb
	fi
	
	for a in pdb/${i}/${i}*.pdb
	do
	    for b in  pdb/${i}/${m}*.pdb
	    do
		A=`basename $a`
		B=`basename $b`
		python3 ~/git/DockQ/DockQ-mod.py -short -useCA  $a $b > pdb/${i}/${A}-${B}.DockQall
	    done
	done
	python3 ~/git/DockQ/DockQ-mod.py -short -useCA pdb/${i}/${i}.pdb pdb/${i}/${m}.pdb > pdb/${i}/${i}.DockQ 
	python3 ~/git/DockQ/DockQ-mod.py -short -useCA pdb/${i}/${i}_reorder.pdb pdb/${i}/${m}.pdb > pdb/${i}/${i}.DockQ.reorder
	cp pdb/${i}/${i}.DockQ pdb/${i}/${i}-${m}.DockQall
	cp pdb/${i}/${i}.DockQ.reorder pdb/${i}/${i}-${m}_reorder.DockQall
	grep -H DockQ pdb/${i}/${i}*.DockQall pdb/${i}/${i}.DockQ pdb/${i}/${i}.DockQ.reorder | sort -nk2 |tail -1 >  pdb/${i}/${i}.DockQall2
	
    
    
	# Original and different chain order
	

	~/git/bioinfo-toolbox/trRosetta/MMalign pdb/${i}/${i}.pdb pdb/${i}/${m}.pdb > pdb/${i}/${i}.MM
	~/git/bioinfo-toolbox/trRosetta/MMalign pdb/${i}/${i}_reorder.pdb pdb/${i}/${m}.pdb > pdb/${i}/${i}.MM.reorder
    fi
    
done <  map_petras.tsv


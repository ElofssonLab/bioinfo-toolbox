#!/bin/bash -x

for i in `cat ../data/pfam/pfam_ids_orthologs_5_single.list`
do
    for j in ../results/uniprot_pfam_annotations/*all.csv
    do
	k=`basename $k all.csv`
	grep $i $j > ../results/uniprot_pfam_annotations/$k.csv
    done
done

#!/bin/bash


echo -n "Eukaryota: " > $1.taxonomy 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Eukaryota >> $1.taxonomy
echo -n "Bacteria: " >> $1.taxonomy 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Bacteria >> $1.taxonomy
echo -n "Archaea: " >> $1.taxonomy 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Archaea >> $1.taxonomy
echo -n "Virus: " >> $1.taxonomy 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Virus >> $1.taxonomy


#!/bin/bash


echo -n "Eukaryota: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Eukaryota 
echo -n "Bacteria: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Bacteria 
echo -n "Archaea: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Archaea 
echo -n "Virus: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab
done | gawk -F"\t" '{print $10}' | grep -c Virus 


#!/bin/bash 


set -euo pipefail
#set +o pipefail
#set +e
#set +u

tempfile=/tmp/kingdoms.$$

echo -n "Eukaryota: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
     grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab || true;  
done  | { gawk -F"\t" '{print $10}' || true; } | { grep -c Eukaryota  || true; }
echo -n "Bacteria: " 

for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab  || true  
done | { gawk -F"\t" '{print $10}'  || true; }  | { grep -c Bacteria  || true  ; }

echo -n "Archaea: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab  || true  
done | { gawk -F"\t" '{print $10}'  || true; } | { grep -c Archaea  || true  ; }

echo -n "Virus: " 
for i in `grep \> $1 | gawk '{print $2}' `
do
    grep -m 1 $i ~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab  || true  
done | { gawk -F"\t" '{print $10}'  || true; }  | { grep -c Virus   || true ; }


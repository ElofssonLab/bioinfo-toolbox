#!/bin/bash

echo "ALL INTERACTIONS"

find tm-pairs/ -name "*\**" -exec rm {} \;

#bin/extractinteractionarea.py  | sort  > interacting-pairs.tsv

echo -n "All pairs: " ; wc -l interacting-pairs.tsv
echo -n "90% pairs: " ; 
for i in `grep \> all_90.seq | sed s/\>//g  | sed s/tm.//g | sed s/.pdb//g `; do grep $i interacting-pairs.tsv ;  done | sort -u > interacting-pairs-90.tsv
wc -l interacting-pairs-90.tsv
echo -n "40% pairs: " ; 
for i in `grep \> all_40.seq | sed s/\>//g  | sed s/tm.//g | sed s/.pdb//g `; do grep $i interacting-pairs.tsv ;  done | sort -u > interacting-pairs-40.tsv
wc -l interacting-pairs-40.tsv


#cdhit -n 2 -c 0.4  -bak 1 -i interacting-pairs-40.seq -o interacting-pairs-40-cdhit.seq
#cdhit  -c 0.9  -bak 1 -i interacting-pairs-90.seq -o interacting-pairs-90-cdhit.seq

bin/filterpairlist.py interacting-pairs-40.tsv interacting-pairs-40-cdhit.seq  | grep _ > interacting-pairs-40-cdhit.tsv
bin/filterpairlist.py interacting-pairs-90.tsv interacting-pairs-90-cdhit.seq  | grep _ > interacting-pairs-90-cdhit.tsv

echo -n "non-identical pairs 90: "
wc -l interacting-pairs-90-cdhit.tsv

echo -n "non-identical pairs 40: "
wc -l interacting-pairs-40-cdhit.tsv




echo "Interactions in TM region"

#bin/extractinteractionareaTM.py  | sort  > interacting-pairs-TM.tsv

gawk '{if ($6>500 && ($4 > 500. && $5 > 500.)) print $0}'  interacting-pairs-TM.tsv | sort > interacting-pairs-TMarea.tsv
gawk '{if ($10>500 && ($4 > 500. && $5 > 500.)) print $0}'  interacting-pairs-TM.tsv | sort  > interacting-pairs-nonTMarea.tsv
gawk '{if ($6>500 && $10>500 && ($4 > 500. && $5 > 500.)) print $0}'  interacting-pairs-TM.tsv | sort  > interacting-pairs-BOTH.tsv


echo -n "All pairs: " ; wc -l interacting-pairs-TMarea.tsv
echo -n "90% pairs: " ; 
for i in `grep \> all_90.seq | sed s/\>//g  | sed s/tm.//g |sed s/\_.//g | sed s/.pdb//g `; do grep $i interacting-pairs-TMarea.tsv ;  done | sort -u > interacting-pairs-TMarea-90.tsv
wc -l interacting-pairs-TMarea-90.tsv
echo -n "40% pairs: " ; 
for i in `grep \> all_40.seq | sed s/\>//g  | sed s/tm.//g | sed s/\_.//g | sed s/.pdb//g `; do grep $i interacting-pairs-TMarea.tsv ;  done | sort -u > interacting-pairs-TMarea-40.tsv
wc -l interacting-pairs-TMarea-40.tsv


#cdhit -n 2 -c 0.4  -bak 1 -i interacting-pairs-40.seq -o interacting-pairs-40-cdhit.seq
#cdhit  -c 0.9  -bak 1 -i interacting-pairs-90.seq -o interacting-pairs-90-cdhit.seq

bin/filterpairlist.py interacting-pairs-TMarea-40.tsv interacting-pairs-40-cdhit.seq  | grep _ > interacting-pairs-TMarea-40-cdhit.tsv
bin/filterpairlist.py interacting-pairs-TMarea-90.tsv interacting-pairs-90-cdhit.seq  | grep _ > interacting-pairs-TMarea-90-cdhit.tsv

echo -n "non-identical pairs 90: "
wc -l interacting-pairs-TMarea-90-cdhit.tsv

echo -n "non-identical pairs 40: "
wc -l interacting-pairs-TMarea-40-cdhit.tsv




# nonTM
echo "Interactions in nonTM region"

echo -n "All pairs: " ; wc -l interacting-pairs-nonTMarea.tsv
echo -n "90% pairs: " ; 
for i in `grep \> all_90.seq | sed s/\>//g  | sed s/tm.//g | sed s/.pdb//g `; do grep $i interacting-pairs-nonTMarea.tsv ;  done | sort -u > interacting-pairs-nonTMarea-90.tsv
wc -l interacting-pairs-nonTMarea-90.tsv
echo -n "40% pairs: " ; 
for i in `grep \> all_40.seq | sed s/\>//g  | sed s/tm.//g | sed s/.pdb//g `; do grep $i interacting-pairs-nonTMarea.tsv ;  done | sort -u > interacting-pairs-nonTMarea-40.tsv
wc -l interacting-pairs-nonTMarea-40.tsv


#cdhit -n 2 -c 0.4  -bak 1 -i interacting-pairs-40.seq -o interacting-pairs-40-cdhit.seq
#cdhit  -c 0.9  -bak 1 -i interacting-pairs-90.seq -o interacting-pairs-90-cdhit.seq

bin/filterpairlist.py interacting-pairs-nonTMarea-40.tsv interacting-pairs-40-cdhit.seq  | grep _ > interacting-pairs-nonTMarea-40-cdhit.tsv
bin/filterpairlist.py interacting-pairs-nonTMarea-90.tsv interacting-pairs-90-cdhit.seq  | grep _ > interacting-pairs-nonTMarea-90-cdhit.tsv

echo -n "non-identical pairs 90: "
wc -l interacting-pairs-nonTMarea-90-cdhit.tsv

echo -n "non-identical pairs 40: "
wc -l interacting-pairs-nonTMarea-40-cdhit.tsv

   


# BOTH
echo "Interactions in BOTH region"

echo -n "All pairs: " ; wc -l interacting-pairs-BOTH.tsv
echo -n "90% pairs: " ; 
for i in `grep \> all_90.seq | sed s/\>//g  | sed s/tm.//g | sed s/.pdb//g `; do grep $i interacting-pairs-BOTH.tsv ;  done | sort -u > interacting-pairs-BOTH-90.tsv
wc -l interacting-pairs-BOTH-90.tsv
echo -n "40% pairs: " ; 
for i in `grep \> all_40.seq | sed s/\>//g  | sed s/tm.//g | sed s/.pdb//g `; do grep $i interacting-pairs-BOTH.tsv ;  done | sort -u > interacting-pairs-BOTH-40.tsv
wc -l interacting-pairs-BOTH-40.tsv


#cdhit -n 2 -c 0.4  -bak 1 -i interacting-pairs-40.seq -o interacting-pairs-40-cdhit.seq
#cdhit  -c 0.9  -bak 1 -i interacting-pairs-90.seq -o interacting-pairs-90-cdhit.seq

bin/filterpairlist.py interacting-pairs-BOTH-40.tsv interacting-pairs-40-cdhit.seq  | grep _ > interacting-pairs-BOTH-40-cdhit.tsv
bin/filterpairlist.py interacting-pairs-BOTH-90.tsv interacting-pairs-90-cdhit.seq  | grep _ > interacting-pairs-BOTH-90-cdhit.tsv

echo -n "non-identical pairs 90: "
wc -l interacting-pairs-BOTH-90-cdhit.tsv

echo -n "non-identical pairs 40: "
wc -l interacting-pairs-BOTH-40-cdhit.tsv

   

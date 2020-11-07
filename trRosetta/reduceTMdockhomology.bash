#!/bin/bash


code=`echo $1  | sed s/.*output_//g  | cut -c 1-4 `

c="XXXX"
# We only check the first 100 hits to save some time.
for i in `grep -vE "^#" $1 |cut -c 1-4  | head -100 ` 
    do
	u=`bl2seq -p blastp -i seq/${code}_u1_A.fasta -j pdbsequences/$i.fasta | grep Expect | head -1 |  gawk '{print $8}'| sed "s/,//g"`
	v=`bl2seq -p blastp -i seq/${code}_u2_A.fasta -j pdbsequences/$i.fasta |  grep Expect | head -1 | gawk '{print $8}'| sed "s/,//g"`
	c+=`echo  $i,$u,$v |gawk -F "," '{if ($2<0.01 && $3 < 0.01) print "|"$1}'` 
done
#echo $c
grep -vE $c $1

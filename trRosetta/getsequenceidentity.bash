#!/bin/bash 

code=`basename $1 _FSA.txt | sed s/output_//g`
tmp=/tmp/bl2seq.out.$$
top=/tmp/tmpfile.out.$$

# We only do it for top 100 hits
grep -v \# $1 | head -100 | sort -u > $top


r='ldkfasklfalfjakljf'
while read -r line
do
    if [[ ${line:0:1} != "#" ]]
    then
	j=`echo $line | gawk '{print $1}' | cut -c 1-4`
	bl2seq -p blastp -e 1000 -D 1 -i pdbsequences/${code}.fasta -j pdbsequences/${j}.fasta -o $tmp
	evalue=`tail -2 $tmp | head -1 | gawk '{ if ( $11>0 )  print $11  ;  else  print '9999'}'`
	seqid=`tail -2 $tmp |  head -1 | gawk '{ if ( $3 != "hits" )  print $3  ;  else  print '0'}'`
	hit=`echo $code $j $evalue $seqid | gawk '{if ($3<0.01) {print $2} }'`
	if [[ $hit > 0 ]]
	then
	    #echo $hit
	    r=$r"|"$hit
	fi
	rm  -f $tmp
    fi
done < $top

echo $r

grep -vE $r $1  > $1.nohomology
rm -f $top

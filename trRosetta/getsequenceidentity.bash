#!/bin/bash 


code=`basename $1 _FSA.txt | sed s/output_//g`
tmp=/tmp/bl2seq.out.$$
top=/tmp/tmpfile.out.$$

# We only do it for top 100 hits
grep -v \# $1 | head -200 | gawk '{print $1}'  | sort -u > $top


r='ldkfasklfalfjakljf'
while read -r line
do
    if [[ ${line:0:1} != "#" ]]
    then
	j=`echo $line | gawk '{print $1}' | cut -c 1-4`
	bl2seq -p blastp -e 1000 -D 1 -i pdbsequences/${code}.fasta -j pdbsequences/${j}.fasta -o $tmp
	evalueA=` grep -E "_1\|.*\|" $tmp | grep -vE ^\# | sort -rgk +11 | tail -1 | gawk '{  if ( $3 >0  )  print $11  ;  else  print '9999'}'`
	evalueB=` grep -E "_2\|.*\|" $tmp  | grep -vE ^\# | sort -rgk +11 | tail -1 | gawk '{  if ( $3 >0  )  print $11  ;  else  print '9999'}'`

	seqidA=`grep -E "_1\|.*\|"  $tmp  | grep -vE ^\# | sort -gk +3 | tail -1  | gawk '{ if ( $3 > 0 )  print $3  ;  else  print '0'}'`
	seqidA=`grep -E "_2\|.*\|"  $tmp  | grep -vE ^\# | sort -gk +3 | tail -1  | gawk '{ if ( $3 > 0 )  print $3  ;  else  print '0'}'`
	hit=`echo $code $j $evalueA $evalueB $seqidA $seqidB | gawk '{if ($3<0.00001 && $4<0.00001 ) {print $2} }'`
	#echo  $code $j $evalue $seqid $hit
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

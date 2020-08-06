#!/bin/bash 

code=`basename $1 _FSA.txt | sed s/output_//g`

echo ${code} , $1
for i in 1 2 3 4 5 6 7 8 9 10
do
    p1=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $26}'`
    p2=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $27}'`
    p3=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $28}'`
    p4=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $29}'`
    p5=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $30}'`
    p6=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $31}'`
    p7=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $32}'`
    p8=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $33}'`
    p9=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $34}'`
    p10=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $35}'`
    p11=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $36}'`
    p12=`grep -vE ^\# $1 | head -$i  | tail -1 | gawk '{print   $37}'`
    #echo $p1,$p2,$p3,$p4,$p5,$p6,$p7,$p8,$p9,$p10,$p11,$p12
    python3 ~/git/bioinfo-toolbox/trRosetta/translatepdb.py -p pdb/${code}_u2_A.pdb -o rotatedpdb/${code}_u2_A_$i.pdb -r " $p1 "  " $p2 "  " $p3 "  " $p4 "  " $p5 "  " $p6 "  " $p7 "  " $p8 "  " $p9 "  " $p10 "  " $p11 "  " $p12 " 
done


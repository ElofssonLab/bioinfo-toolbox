#!/bin/bash -x

find ./pdb/ -size 0 -exec rm {} \;
for i in `find pdb/ -name "*DockQ" -exec grep -Hc already {} \; | grep -v :0 |sed s/:.*//g`
do
    rm $i
done

for i in pdb/*/*DockQ
do
    j=`basename $i .DockQ`
    echo -n $j ","
    gawk '{print $2}' $i
done | sed "s/ //g"  > DockQ.csv

for i in `find pdb/ -name "*DockQ.reorder" -exec grep -Hc already {} \; | grep -v :0 |sed s/:.*//g`
do
    rm $i
done

for i in pdb/*/*DockQ.reorder
do
    j=`basename $i .DockQ.reorder`
    echo -n $j ","
    gawk '{print $2}' $i
done | sed "s/ //g"  > DockQreorder.csv

for i in pdb/*/*DockQall2
do
    j=`basename $i .DockQall2`
    echo -n $j ","
    gawk '{print $2}' $i
done | sed "s/ //g"  > DockQall.csv


# A few MM files crashed (as they are only CA)
find pdb/ -name "*.MM" -size 1 -exec rm {} \;

for i in pdb/*/*MM
do
    j=`basename $i .MM`
    echo -n $j ","
    grep TM-score= $i | head -2 | sort -n -k 2 | tail -1 |gawk '{print $2}'
done  |sed "s/ //g" >  MM.csv

find pdb/ -name "*.MM.reorder" -size 1 -exec rm {} \;
for i in pdb/*/*MM.reorder
do
    j=`basename $i .MM`
    echo -n $j ","
    grep TM-score= $i | head -2 | sort -n -k 2 | tail -1 |gawk '{print $2}'
done  |sed "s/ //g" >  MMreorder.csv

find pdb/ -name "*.MMall" -size 1 -exec rm {} \;
find pdb/ -name "*.MMall2" -size 0 -exec rm {} \;
for i in pdb/*/*MMall2
do
    j=`basename $i .MMall2 `
    a=`grep TM-score= $i | sed "s/.*= //g" | sed "s/\ .*//g" | sort -n | tail -1 `
    b=`grep TM-score= pdb/$j/$j.MM* | head -2 | sort -n -k 2 | tail -1 |gawk '{print $2}'`
    echo -n $j ","
    echo $a $b | gawk '{ if($1>$2){print $1}else{print $2}}'
done  |sed "s/ //g" > MMall.csv

#for i in pdb/*/*MMall2
#do
#    #j=`basename $i .MMall | sed s/\-[0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9]//g`
#    j=`basename $i .MMall2 `
#    a=`grep TM-score= $i | head -2 | sort -n -k 5 | tail -1 |sed "s/.*= //g" | sed "s/ .*//g' `
#    b=`grep TM-score= pdb/$j/$j.MM* | head -2 | sort -n -k 2 | tail -1 |gawk '{print $2}'`
#    echo -n $j ","
#    echo $a $b | gawk '{ if($1>$2){print $1}else{print $3}}'
#    #grep TM-score= $i | head -2 | sort -n -k 5 | tail -1 |gawk '{print $2","$3","$4","$5}' 
#done  |sed "s/ //g" >  MMall.csv

for i in pdb/*/*.pLDDT
do
    k=`dirname $i`
    j=`basename $k`
    echo -n $j ","
    #gawk '{if ($1=="IF_NumRes:") {N=$2}  else if ($1=="IF_pLDDT") {I=$2}   else {print N ","  I "," $2}}' $k/$j.pLDDT
    #grep Summary: $k/$j.pLDDT | sed "s/\[//g" | sed "s/\]//g" | sed "s/Summary://g"
    gawk '{if ($1=="IF_NumRes:") {N=$2}  
    	 else if ($1=="IF_pLDDT") {I=$2}   
    	 else if ($1=="pLDDT") {p=$2}   
    	 else if ($1=="Summary:") {print N ","  I ","  p  }}'    $k/$j.pLDDT 
done |sed "s/ //g" > pLDDT.csv

echo "Name,NumRes,IF_plDDT,plDDT" > negatome-pLDDT.csv
for i in negatome-pLDDT/*
do
    j=`basename $i .pLDDT`
    echo -n $j ","
    #gawk '{if ($1=="IF_NumRes:") {N=$2}  else if ($1=="IF_pLDDT") {I=$2}   else {print N ","  I "," $2}}' $k/$j.pLDDT
    #grep Summary: $k/$j.pLDDT | sed "s/\[//g" | sed "s/\]//g" | sed "s/Summary://g"
    gawk '{if ($1=="IF_NumRes:") {N=$2}  
    	 else if ($1=="IF_pLDDT") {I=$2}   
    	 else if ($1=="pLDDT") {p=$2}   
    	 else if ($1=="Summary:") {print N ","  I ","  p }}'    $i
done |sed "s/ //g" >> negatome-pLDDT.csv

bin/mergecsv.py

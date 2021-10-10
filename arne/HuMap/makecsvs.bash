#!/bin/bash -x

find ./pdb/ -size 0 -exec rm {} \;

for i in pdb/*/*DockQ
do
    j=`basename $i .DockQ`
    echo -n $j ","
    gawk '{print $2}' $i
done | sed "s/ //g"  > DockQ.csv


# A few MM files crashed (as they are only CA)
find pdb/ -name "*.MM" -size 1 -exec rm {} \;

for i in pdb/*/*MM
do
    j=`basename $i .MM`
    echo -n $j ","
    grep TM- $i | head -2 | sort -n -k 2 | tail -1 |gawk '{print $2}'
done  |sed "s/ //g" >  MM.csv

for i in pdb/*/unrelaxed_model_1.pdb
do
    k=`dirname $i`
    j=`basename $k`
    echo -n $j ","
    #gawk '{if ($1=="IF_NumRes:") {N=$2}  else if ($1=="IF_pLDDT") {I=$2}   else {print N ","  I "," $2}}' $k/$j.pLDDT
    #grep Summary: $k/$j.pLDDT | sed "s/\[//g" | sed "s/\]//g" | sed "s/Summary://g"
    gawk '{if ($1=="IF_NumRes:") {N=$2}  
    	 else if ($1=="IF_pLDDT") {I=$2}   
    	 else if ($1=="pLDDT") {p=$2}   
    	 else if ($1=="if_plddt_sum") {s=$2}   
    	 else if ($1=="Summary:") {print N ","  I ","  p ","  s  }}'    $k/$j.pLDDT 
done |sed "s/ //g" > pLDDT.csv

bin/mergecsv.py

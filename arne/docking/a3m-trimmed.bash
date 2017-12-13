#!/bin/bash

for i in *seq/*a3m
do
    j=`basename $i .a3m`
    d=`dirname $i`
    c=`echo $j | sed s/.fa.HH.*//g `
    ~/git/PconsC3/extra/arne/MSA/a3mToTrimmed.py -o -name $c $i > $d/$j.trimmed
done

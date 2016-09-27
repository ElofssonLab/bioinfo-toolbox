#!/bin/bash

INDIR=$1

# obtained from: http://stephan.paukner.cc/syslog/archives/362-Averaging-an-image-sequence-with-ImageMagick.html

i=0
for file in ${INDIR}/*.png; do
    echo -n "$file.. "
    if [ $i -eq 0 ]; then
        cp $file ${INDIR}/avg.png
    else
        convert $file ${INDIR}/avg.png -fx "(u+$i*v)/$[$i+1]" ${INDIR}/avg.png
    fi
    i=$[$i+1]
done

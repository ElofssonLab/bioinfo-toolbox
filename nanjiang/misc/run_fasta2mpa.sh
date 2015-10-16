#!/bin/bash

# Description:
# convert mfa file to mpa file and then clean mfa file
errfile=$(mktemp /tmp/tmp.run_fasta2mpa.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  

trap 'rm -f "$errfile"' INT TERM EXIT

while [ "$1" != "" ]; do
    infile="$1"
    stemfilename=${infile%.*}
    $DATADIR3/wk/MPTopo/src/msa2mpa.py $infile > $stemfilename.mpa 2> $errfile
    if [ ! -s $errfile -a -s "$stemfilename.mpa" ]; then
        echo "$infile ==> $stemfilename.mpa"
        rm -f $infile
    else
        cat $errfile  >&2
    fi
    shift
done

rm -f $errfile

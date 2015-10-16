#!/bin/bash

# check ok on 2013-04-18
# nanjiang@casio /data3/wk/MPTopo/pfamAna_refpro/cellular_filter_fragment/msa_kalignp_splitclanseq_filterNonTM
# $ /data3/wk/MPTopo/src/tmp_checkmpa.sh `find *.msa.fa | randlist | head -n 10`

binpath=$DATADIR3/wk/MPTopo/src

CheckMPA() {
    local infile=$1
    local mpafile=$infile.mpa
    $binpath/msa2mpa.py $infile -o $mpafile

    local rewritten_mfafile=$infile.rewritten
    $binpath/rewritefasta.py $infile -o $rewritten_mfafile

    local converted_mfafile=$infile.convertedfrommpa
    $binpath/mpa2msa.py $mpafile -o $converted_mfafile

    msg=`diff $rewritten_mfafile $converted_mfafile | head -n 1`
    if [ "$msg" != "" ]; then
        echo failed for $infile
    else
        echo OK for $infile
    fi
}
while [ "$1" != "" ]; do
    CheckMPA "$1"
    shift
done

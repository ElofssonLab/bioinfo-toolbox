#!/bin/bash

rundir=`dirname $0`

cd $rundir

# Step 1. backup source codes
bash ./bakfile.sh

# Step 2, install everything to ../bin
OPT1="--exclude=~* --exclude=bak --exclude=*~ --exclude=.*.sw[mopn] --exclude=*conflicted*copy "
rsync -auvz $OPT1 --exclude="install.sh" ./ ../bin/

# Step 3, install to binpath 
scriptlist="
myfunc.py
mybase.py
libtopologycmp.py
getfastaid.py
catfasta.py
randfasta.py
mydb_common.py
my_formatdb.py
my_extractdb.py
my_indexformatconvert.py
splitfasta.py
countseq.py
selectLineByID.py
msa2mpa.py
mpa2msa.py
getseqlen.py
splitfolder.py
indexfasta.py
"

for f in $scriptlist; do
    /bin/rm -f $BINPATH/$f
    /bin/cp -f $f $BINPATH/
done

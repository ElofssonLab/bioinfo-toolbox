#!/bin/bash -x

# Just a script that actually checks if everything is nicely packes in the way it should be



# Then we need to delete all directories that we already have as tar files.

dir=`pwd`

#for i in `find ./PF00001.18/ -type d -maxdepth 1 -mindepth 1 -type d `
for i in `find ./  -maxdepth 2 -mindepth 2 -type d `
do
    j=`basename $i`
    k=`dirname $i`
    cd $k
    if [ -s $j.tar.gz ]
    then
	tar -zxvf $j.tar.gz
    fi
    tar -zcvf ${j}.tar.gz ${j}  --remove-files  
    cd ${dir}
done

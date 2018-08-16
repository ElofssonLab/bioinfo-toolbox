#!/bin/bash -lx
#SBATCH --output=parse.%A_%a.out
#SBATCH --error=parse.%A_%a.out
#SBATCH --array=1-335
#SBATCH -c 1
#SBATCH -t 30:00
#SBATCH -A SNIC2016-10-22




for id in PF*
do
    dir=`pwd`/$id
    
    scratch=$SNIC_TMP/arnee/Summary/$id/
    mkdir -p $scratch
    cd $scratch
    
    
    sleep 2 # waiting for filesystem
    ln -fs $dir/*.out ./
    ln -fs $dir/*.raw ./
    ln -fs $dir/*.gneff ./
    
    for i in $dir/*cm.tar.gz # $dir/conf*[04].tar.gz # $dir/*_mem.tar.gz
    do
	j=`basename $i .tar.gz`
	if [ ! -s $dir/${j}_summary.csv ]
	then
	    if [ -s $dir/${j}_proq3.tar.gz ]
	    then
		tar -zxf $dir/${j}_proq3.tar.gz
		sleep 1
		$dir/../bin/GetAllScores.py ${j}.tar.gz >  $dir/${j}_summary.csv
		#	mv ${j}*out $dir/
		sleep 1
		rm -r ${j}/
	    fi
	fi
    done
    
    cd $dir/../
done

    

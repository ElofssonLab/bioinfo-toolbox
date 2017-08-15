#!/bin/bash -x
#SBATCH -c 2
#SBATCH -t 00:00:01
#SBATCH -A SNIC2016-10-22


cat /proc/$$/status | grep Cpus_allowed_list | gawk '{print $2}' | sed "s/-/ /g" | gawk '{if ($2>0){print 1+$2-$1}else{print 1}}'

#!/bin/bash -l
#SBATCH --job-name=filter_vall
#SBATCH --output=filter_vall.slurm.out
#SBATCH --error=filter_vall.slurm.error
#SBATCH -n 1
#SBATCH -t 50:00:00
OUTFILE=$1
DB=$2
ID_lower=$3
python ~/toolbox/rosetta/filter_vall.py ${OUTFILE} ${DB} ${ID_lower}

#!/bin/bash -l
#SBATCH --job-name=python
#SBATCH --output=python.slurm.out
#SBATCH --error=python.slurm.error
#SBATCH -A b2011166
#SBATCH -p node
#SBATCH -N 1 -n 1
#SBATCH -t 10:00:00
INFILE=$1
NAME=$2
OUTFILE=$3
module load bioinfo-tools
python analyse_alignment.py ${INFILE} ${NAME} > ${OUTFILE}

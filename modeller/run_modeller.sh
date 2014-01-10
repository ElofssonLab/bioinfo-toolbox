#!/bin/bash -l
#SBATCH --output=modeller.slurm.out
#SBATCH --error=modeller.slurm.error
#SBATCH -n 1
#SBATCH -t 5:00:00
ID=$1
SEQFILE=$2
TEMPL_ID=$3
TEMPL_PDB=$4
CHAIN=$5
/home/x_mirmi/bin/modeller9.12/bin/modpy.sh python /home/x_mirmi/bioinfo-toolbox/modeller/align_to_template.py ${ID} ${SEQFILE} ${TEMPL_ID} ${TEMPL_PDB} ${CHAIN} ${ID}_${TEMPL_ID}.ali
/home/x_mirmi/bin/modeller9.12/bin/modpy.sh python /home/x_mirmi/bioinfo-toolbox/modeller/model_simple.py ${ID}_${TEMPL_ID}.ali ${ID} ${TEMPL_ID}

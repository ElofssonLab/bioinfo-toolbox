#!/bin/bash -x
#SBATCH -A SNIC2018-1-26
#SBATCH --output=annotate.%A_%a.out
#SBATCH --error=annotate.%A_%a.out
#SBATCH --array=1-50
#SBATCH -J annotate_uniprot.%A_%a
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -c 8 # Or 14, or any other number of cores to use for tblastn
# -c can also be specified as arg to "sbatch" and overrides the value
# specified here.

source /pfs/nobackup/home/w/wbasile/venv-work/bin/activate

python ./annotate_protein_fasta.py

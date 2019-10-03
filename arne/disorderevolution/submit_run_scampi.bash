#!/bin/bash -x
#SBATCH -A SNIC2019-3-319
#SBATCH --output=scampi.%A_%a.out
#SBATCH --error=scampi.%A_%a.out
#SBATCH --array=1-50
#SBATCH -J scampi.%A_%a
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -c 1 # Or 14, or any other number of cores to use for tblastn
# -c can also be specified as arg to "sbatch" and overrides the value
# specified here.


source /pfs/nobackup/home/w/wbasile/venv-work/bin/activate

python /pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/bin/run_iupred.py

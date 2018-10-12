#!/bin/bash
#SBATCH -A SNIC2018-1-26
#SBATCH -J create-annotate-single
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -c 8 # Or 14, or any other number of cores to use for tblastn
# -c can also be specified as arg to "sbatch" and overrides the value
# specified here.

source /pfs/nobackup/home/w/wbasile/venv-work/bin/activate

#python ./do_blast.py
python create_final_dataset_single.py

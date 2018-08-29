#!/bin/bash
#SBATCH -A SNIC2017-11-7
#SBATCH -J print_prop
#SBATCH -t 00:30:00
#SBATCH -n 1
#SBATCH -c 8 # Or 14, or any other number of cores to use for tblastn
#SBATCH --mem=32GB
# -c can also be specified as arg to "sbatch" and overrides the value
# specified here.

source ../../venv-work/bin/activate
python ./print_average_prop_per_protein.py

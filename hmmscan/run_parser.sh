#!/bin/bash -l
#SBATCH --job-name=pfparse
#SBATCH --output=hmmer3.slurm.out
#SBATCH --error=hmmer3.slurm.error
#SBATCH -A b2011166
#SBATCH -p node
#SBATCH -N 1 -n 8
#SBATCH -t 00:45:00
python testscript.py ../databases/Pfam/current_release/Pfam-A.full > dict_acc_numseq.py

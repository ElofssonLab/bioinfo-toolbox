#!/bin/bash -l
#SBATCH --job-name=hmmer3
#SBATCH --output=hmmer3.slurm.out
#SBATCH --error=hmmer3.slurm.error
#SBATCH -A b2011166
#SBATCH -p node
#SBATCH -N 1 -n 8
#SBATCH -t 35:00:00
module load openmpi
/sw/apps/bioinfo/hmmer/3.0_kalkyl_intel/bin/hmmscan --domtblout pdb_seqres_2013-02_2.scan ~/glob/databases/Pfam/current_release/Pfam-A.hmm ~/glob/databases/PDB/current_release/pdb_seqres.txt

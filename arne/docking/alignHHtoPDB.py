#!/usr/bin/env python

use Bio

# Make profile with HHblits/jackhmmer

runhhblits.py seqA.fa
runhbblits.py seqB.fa

#BUild HMM

# Align HMM to PDB

align seqA.hmm pdbDB
align seqA.hmm pdbDB

# Extract PDB files.
for i in align
    if seqA.pdb == seqB.pdb
    extractpdb seqA
    extractpdb seqB
# Check surface areas
    freesasa pdbA
    freesasa pdbB
    freesasa pdbA+B
    #Squirl ??
    if (ab-a-b<cutoff) print pdbA+B
        

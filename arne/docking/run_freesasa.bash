#!/bin/bash
for i in tm-chains/*.pdb
do
    /pfs/nobackup/home/a/arnee/Software/freesasa-2.0.2/src/freesasa --format pdb $i > $i.pdb-Bvalue.pdb
done

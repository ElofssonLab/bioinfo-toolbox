#!/bin/bash -x


# We actually need to renumber residues first
A=/tmp/temp1.pdb.$$
B=/tmp/temp2.pdb.$$

#pdb_reres $1 > $A 
#pdb_reres $2 > $B
#python3 ~/git/DockQ/DockQ-mod.py -short -useCA  $A $B
python3 ~/git/DockQ/DockQ-mod.py -short -useCA  $1 $2

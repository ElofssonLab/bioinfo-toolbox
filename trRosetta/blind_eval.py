from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import argparse

if __name__ == "__main__":
    p = argparse.ArgumentParser(description = '- Merging three NPZ files (each file + interaction area)')
    p.add_argument('-s', required= True, help='complex structure')
    p.add_argument('-c', required= False, default=3, help='clash thr distance in A')
    ns = p.parse_args()

    p = PDBParser(QUIET=True)
    dimer = p.get_structure('', ns.s)

    rA = Selection.unfold_entities(dimer[0]['A'], 'R')
    rB = Selection.unfold_entities(dimer[0]['B'], 'R')
    mind = ''
    clashes = 0
    for r1 in rA:
        for r2 in rB:
            d12 = r1['CA']-r2['CA']
            if d12 < 3: clashes += 1
            if mind == '' or d12 <= mind: mind = d12 

    print ('Complex: {} Clash#: {} Chain_distance: {}'.format(ns.s, clashes, mind))


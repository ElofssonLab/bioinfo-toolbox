from Bio.PDB import *
from Bio.PDB.DSSP import DSSP

if __name__ == "__main__":

    p.add_argument('-s', required= True, help='complex structure')
    p.add_argument('-c', required= False, default=3, help='clash thr distance in A')
    ns = p.parse_args()

    p = PDBParser(QUIET=True)
    dimer = p.get_structure('', ns.s)

    atomsA = Selection.unfold_entities(dimer[0]['A'], 'A')
    atomsB = Selection.unfold_entities(dimer[0]['B'], 'A')
    clashes = [[a, b] for a in atomsA for b in atomsB if a-b < 3]
    mindistance = min([a-b for a in atomsA for b in atomsB])
    print ('Clash#: {} Chain_distance: {}'.format(len(clashes), mindistance))


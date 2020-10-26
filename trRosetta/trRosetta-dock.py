import sys,os,json
import tempfile
import argparse
import numpy as np
import random as rnd
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose, renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.protocols.rigid import *
from foldclass import constr_fold
from dockclass import constr_dock

def main():
    
    trf = constr_fold(ns.src, ns.dat)
    trd = constr_dock(ns.src, ns.dat)

    ##### Pose initialization #####
    p1 = trf.pose_from_fasta(ns.fasta1)
    p2 = trf.pose_from_fasta(ns.fasta2)
    trf.par['LENSEQ1'] = trd.par['LENSEQ1'] = len(p1.sequence())
    trf.par['LENSEQ2'] = len(p2.sequence())
    trf.par['SEQ'] = trd.par['SEQ'] = p1.sequence()+p2.sequence()
    trf.par['LENSEQ'] = trf.par['LENSEQ1']+trf.par['LENSEQ2']
 
    trf.mutate_gly(p1.sequence(), p1, 'ALA')
    trf.mutate_gly(p2.sequence(), p2, 'ALA')

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(True)

    min_mover_vdw = MinMover(mmap, trf.sfs[2], 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover_vdw.max_iter(500)
    trf.remove_clashes(p1, min_mover_vdw, trf.sfs[2])

    npz = np.load(ns.npz)
    if npz['dist'].shape[0] == trf.par['LENSEQ']+20:
        sliced = {}
        for key in npz.files:
            key = str(key)
            sliced[key] = np.delete(npz[key], slice(trf.par['LENSEQ'], trf.par['LENSEQ']+21), axis=0)
            sliced[key] = np.delete(sliced[key], slice(trf.par['LENSEQ'], trf.par['LENSEQ']+21), axis=1)
        npz = sliced
    print (npz['dist'].shape, trf.par['LENSEQ'])

    npz1 = trf.npz_selection(npz, 1)
    npz2 = trf.npz_selection(npz, 2)
    rst1 = trf.format_rst(npz1)
    rst2 = trf.format_rst(npz2)

    trf.folding(p1, rst1, mmap)
    array = trf.select_rst(rst1, 1, trf.par['LENSEQ1'], p1.sequence())
    trf.relax(p1, array)

    trf.folding(p2, rst2, mmap)
    array = trf.select_rst(rst2, 1, trf.par['LENSEQ2'], p2.sequence())
    trf.relax(p2, array)

    append_pose_to_pose(p1, p2)
    renumber_pdbinfo_based_on_conf_chains(p1)
    trf.mutate_gly(trf.par['SEQ'], p1, 'GLY')

    #rst = trf.format_rst(npz)
    #array = trf.select_rst(rst, 1, trf.par['LENSEQ'], trf.par['SEQ'], inter=False)
    #array = trd.add_site(npz)
    #trf.relax(p1, array, bb=False)

    p1.dump_pdb(ns.out)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = "Class to implement trRosetta")
    parser.add_argument("npz", type=str, help="input distograms and anglegrams (NN predictions)")
    parser.add_argument("fasta1", type=str, help="input sequence 1")
    parser.add_argument("fasta2", type=str, help="input sequence 2")
    parser.add_argument("out", type=str, help="output model (in PDB format)")
    parser.add_argument("-dat", required=True, type=str, help="folder containing data subfolders (seq & npz)")
    parser.add_argument("-src", required=True, type=str, help="folder containing accessory files")
    ns = parser.parse_args()

    main()



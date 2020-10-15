import sys,os,json
import tempfile
import argparse
import numpy as np
import random as rnd
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose, renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.protocols.rigid import *

os.environ["OPENBLAS_NUM_THREADS"] = "1"

class trRosetta_fold():

    def __init__(self, scriptdir, datadir):
        self.scriptdir = scriptdir.rstrip('/')+'/'
        self.datadir = datadir.rstrip('/')+'/'
        self.tmpdir = tempfile.TemporaryDirectory(prefix=datadir)

        with open(scriptdir+'/data/params.json') as jsonfile:
            self.par = json.load(jsonfile)

        init('-hb_cen_soft -relax:default_repeats 5 \
              -default_max_cycles 200 -out:level 100')

    def pose_from_fasta(self, fasta_file):
        self.sequence = ''
        for line in open(self.datadir+fasta_file): 
            if not line.startswith('>'): self.sequence += line.rstrip()
        pose = pose_from_sequence(self.sequence, 'centroid')
        return pose

    def set_pose_info(self, pose, start, end, chain):
        count = 0
        sequence = pose.sequence()
        for i,a in enumerate(sequence): 
            if i < start: continue
            if i == end: break
            print (i, chain)
            pose.pdb_info().chain(i+1, chain)

    def remove_clashes(self, pose, mover, rc_sf, thr=10, repeat=5):
        for _ in range(0, repeat):
            if float(rc_sf(pose)) < thr: break
            mover.apply(pose)

    def mutate_gly(self, sequence, pose, mut):
        for i,a in enumerate(sequence):
            if a == 'G':
                mutator = rosetta.protocols.simple_moves.MutateResidue(i+1, mut.upper())
                mutator.apply(pose)

    def ISPRED_site(self, ispred, res, thr=0.1):
        site = []
        non_site = []
        for line in open(ispred):
            res += 1
            score = float(line.split()[-1])
            if self.par['SEQ'][res-1] == 'G': continue
            if score > thr: site.append(res)
            else: non_site.append(res)
        return site, non_site

    def ISPRED_top_site(self, ispred, res, top=3):
        site = []
        non_site = []
        for line in open(ispred):
            res += 1
            score = float(line.split()[-1])
            if self.par['SEQ'][res-1] == 'G': continue
            if float(score) == 0.0: non_site.append(res)
            else: 
                for p, r in enumerate(site):
                    if score > r[1]: 
                        site = site[:p]+[[res,score]]+site[p:]
                        break
                if [res,score] not in site: site.append([res,score])           
        site = [el[0] for el in site][:top]

        return site, non_site

    def format_ISPRED_rst(self, siteA, siteB, repA, repB):
        array = ['AtomPair CB {} CB {} FLAT_HARMONIC 8 20 4'.format(a, b) for a in siteA for b in siteB]
        array += ['AtomPair CB {} CB {} FLAT_HARMONIC 1000 10 980'.format(a, b) for a in repA for b in repB]
        array += ['AtomPair CB {} CB {} FLAT_HARMONIC 1000 10 990'.format(a, b) for a in repA for b in siteB]
        array += ['AtomPair CB {} CB {} FLAT_HARMONIC 1000 10 990'.format(a, b) for a in repB for b in siteA]
        return array

    def apply_rst(self, pose, array, constraints_file):
        with open(constraints_file,'w') as f:
            for line in array: f.write(line+'\n')
        constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
        constraints.constraint_file(constraints_file)
        constraints.add_constraints(True)
        constraints.apply(pose)
        os.remove(constraints_file)

    def custom_docking(self, pose, isfile1, isfile2):
        print (ns.out+'_1.pdb')

        mmap = MoveMap()
        mmap.set_bb(False)
        mmap.set_chi(False)
        mmap.set_jump(True)

        sf_fa = create_score_function('docking')
        sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)

        relax = rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sf_fa)
        relax.dualspace(True)
        relax.set_movemap(mmap)

        rotation = 90
        translation = 10
        dock_pert = RigidBodyPerturbMover(1, translation, rotation)

        to_fullatom = SwitchResidueTypeSetMover('fa_standard')
        to_fullatom.apply(pose)

        siteA, repA = self.ISPRED_top_site(isfile1, 0)
        siteB, repB = self.ISPRED_top_site(isfile2, self.par['LENSEQ1'])
        array = self.format_ISPRED_rst(siteA, siteB, repA, repB)
        print ('Extracted {} constraints'.format(len(array)))

        #real = pose_from_pdb('./3f1p/3f1p.pdb')
        #self.apply_rst(real, array, self.tmpdir.name+'/minimize.cst')

        #decoy = Pose()
        for n in range(5):
            #decoy.assign(pose)
            pose.remove_constraints()
            dock_pert.apply(pose)
            self.apply_rst(pose, array, self.tmpdir.name+'/minimize.cst')
            print ('Pose energy before dock:'+str(sf_fa(pose)))
            #pose.dump_pdb(ns.out+'_'+str(n+1)+'_p.pdb')
            relax.apply(pose)
            print (sf_fa(pose))
            pose.dump_pdb(ns.out+'_'+str(n+1)+'.pdb')

    def standard_docking(self, pose, isfile1, isfile2): 
        jobs=1
        dock_jump = 1
        rotation = 15
        translation = 3
        job_output = 'dock_output'
        
        to_centroid = SwitchResidueTypeSetMover('centroid')
        to_fullatom = SwitchResidueTypeSetMover('fa_standard')

        protocols.docking.setup_foldtree(pose, 'A_B', Vector1([dock_jump]))

        sf_low = create_score_function('interchain_cen')
        sf_low.set_weight(rosetta.core.scoring.site_constraint, 10)
        sf_high = create_score_function('docking')
        sf_high.set_weight(rosetta.core.scoring.site_constraint, 10)
        sf_high_min = create_score_function('docking','docking_min')
        sf_high_min.set_weight(rosetta.core.scoring.site_constraint, 10)

        randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump, partner_upstream)
        randomize_downstream = RigidBodyRandomizeMover(pose, dock_jump, partner_downstream)
        dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
        spin = RigidBodySpinMover(dock_jump)
        slide_into_contact = protocols.docking.DockingSlideIntoContact(dock_jump)
        
        movemap = MoveMap()
        movemap.set_jump(dock_jump, True)

        minmover = MinMover()
        minmover.movemap(movemap)
        minmover.score_function(sf_high_min)


        dock_prot = protocols.docking.DockingProtocol()
        dock_prot.set_movable_jumps(Vector1([dock_jump]))
        dock_prot.set_lowres_scorefxn(sf_low)
        dock_prot.set_highres_scorefxn(sf_high_min)

        #docking_low = protocols.docking.DockingLowRes()
        #docking_low.set_movable_jumps(Vector1([dock_jump]))
        #docking_low.set_scorefxn(sf_low)

        #docking_high = protocols.docking.DockingHighRes()
        #docking_high.set_movable_jumps(Vector1([dock_jump]))
        #docking_high.set_scorefxn(sf_high)


def main():
    
    trf = trRosetta_fold(ns.s, ns.d)

    ##### Pose initialization #####
    p1 = pose_from_pdb(ns.str1)
    p2 = pose_from_pdb(ns.str2)
    trf.par['LENSEQ1'] = len(p1.sequence())
    trf.par['LENSEQ2'] = len(p2.sequence())
    append_pose_to_pose(p1,p2)
    renumber_pdbinfo_based_on_conf_chains(p1)
    trf.par['SEQ'] = p1.sequence()
    trf.par['LENSEQ'] = len(p1.sequence())

    print (p1.pdb_info())
    ##### Relaxation #####
    print ('Relaxing ...')

    trf.custom_docking(p1, ns.i1, ns.i2)

    #p1.dump_pdb(ns.out)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = "Class to implement trRosetta")
    parser.add_argument("-npz", type=str, help="input distograms and anglegrams (NN predictions)")
    parser.add_argument("str1", type=str, help="input structure 1")
    parser.add_argument("str2", type=str, help="input structure 2")
    parser.add_argument("out", type=str, help="output model (in PDB format)")
    parser.add_argument("-i1", type=str, help="ISPRED file for prot 1")
    parser.add_argument("-i2", type=str, help="ISPRED file for prot 2")
    parser.add_argument("-s", type=str, help="script directory")
    parser.add_argument("-d", type=str, help="data directory")
    ns = parser.parse_args()
    main()



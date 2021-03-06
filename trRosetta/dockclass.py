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

class constr_dock():

    def __init__(self, scriptdir, datadir):
        self.scriptdir = scriptdir.rstrip('/')+'/'
        self.datadir = datadir.rstrip('/')+'/'
        self.tmpdir = tempfile.TemporaryDirectory(prefix=datadir)

        with open(scriptdir+'/data/params.json') as jsonfile:
            self.par = json.load(jsonfile)

        init('-hb_cen_soft -relax:default_repeats 5 \
              -default_max_cycles 200 -out:level 100')

    def add_site(self, npz):
        array = []
        dist = npz['dist'][:,:,5:]
        sep = self.par['LENSEQ1']
        pcut = self.par['PCUT']
        psumup = np.sum(dist, axis=-1)
        i,j = np.where(psumup>pcut+0.5)
        for a,b in zip(i,j):
            if a>b or abs(a-b)<=5 or np.sign(a-sep)==np.sign(b-sep): continue
            if self.par['SEQ'][a] == 'G' or self.par['SEQ'][b] == 'G': continue
            flatcenter = (np.where(dist[a,b]==np.amax(dist[a,b]))[0][0]//2)+4.5
            rst_line = 'AtomPair CB {a} CB {b} FLAT_HARMONIC {f} 2 1'.format(a=a+1, b=b+1, f=flatcenter)
            if rst_line not in array: array.append(rst_line)
        return array

    def ISPRED_site(self, ispred, res, thr=0.3):
        site = []
        non_site = []
        for line in open(ispred):
            res += 1
            score = float(line.split()[-1])
            if self.par['SEQ'][res-1] == 'G': continue
            if score > thr: site.append(res)
            else: non_site.append(res)
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

        siteA, repA = self.ISPRED_site(isfile1, 0)
        siteB, repB = self.ISPRED_site(isfile2, self.par['LENSEQ1'])
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
            print ('To dock:{}, Real:{}'.format(sf_fa(pose), sf_fa(real)))
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



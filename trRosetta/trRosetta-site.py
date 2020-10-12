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

    def format_rst(self, npz_file):
        #### note about constraint functions: when Y gets small (or negative) it seems to be rewarding.
        #### Example at: https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/constraint-file
        #### look at "Sample Files" section.
        array = []
        SEQ = self.par['SEQ']
        NRES = self.par['LENSEQ']
        DCUT = self.par['DCUT']
        PCUT = self.par['PCUT']
        ASTEP = np.deg2rad(self.par['ASTEP'])
        DSTEP = self.par['DSTEP']
        ALPHA = self.par['ALPHA']
        EBASE = self.par['EBASE']
        EREP = self.par['EREP']
        DREP = self.par['DREP']
        MEFF = self.par['MEFF']

        npz = np.load(self.datadir+npz_file)
        dist,omega,theta,phi = npz['dist'],npz['omega'],npz['theta'],npz['phi']
        rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : []}
        ########################################################
        # dist: 0..20A
        ########################################################
        nres = dist.shape[0]
        bins = np.array([4.25+DSTEP*i for i in range(32)])
        prob = np.sum(dist[:,:,5:], axis=-1)
        bkgr = np.array((bins/DCUT)**ALPHA)
        attr = -np.log((dist[:,:,5:]+MEFF)/(dist[:,:,-1][:,:,None]*bkgr[None,None,:]))+EBASE
        repul = np.maximum(attr[:,:,0],np.zeros((nres,nres)))[:,:,None]+np.array(EREP)[None,None,:]
        dist = np.concatenate([repul,attr], axis=-1)
        bins = np.concatenate([DREP,bins])
        i,j = np.where(prob>PCUT)
        prob = prob[i,j]
        nbins = 35
        step = 0.5
        for a,b,p in zip(i,j,prob):
            if b>a:
                name=self.tmpdir.name+"/%d.%d.txt"%(a+1,b+1)
                with open(name, "w") as f:
                    f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                    f.write('y_axis'+'\t%.3f'*nbins%tuple(dist[a,b])+'\n')
                    f.close()
                rst_line = 'AtomPair %s %d %s %d SPLINE TAG %s 1.0 %.3f %.5f'%('CB',a+1,'CB',b+1,name,1.0,step)
                rst['dist'].append([a,b,p,rst_line])
        print("dist restraints:  %d"%(len(rst['dist'])))
        ########################################################
        # omega: -pi..pi
        ########################################################
        nbins = omega.shape[2]-1+4
        bins = np.linspace(-np.pi-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
        prob = np.sum(omega[:,:,1:], axis=-1)
        i,j = np.where(prob>PCUT)
        prob = prob[i,j]
        omega = -np.log((omega+MEFF)/(omega[:,:,-1]+MEFF)[:,:,None])
        omega = np.concatenate([omega[:,:,-2:],omega[:,:,1:],omega[:,:,1:3]],axis=-1)
        for a,b,p in zip(i,j,prob):
            if b>a:
                name=self.tmpdir.name+"/%d.%d_omega.txt"%(a+1,b+1)
                with open(name, "w") as f:
                    f.write('x_axis'+'\t%.5f'*nbins%tuple(bins)+'\n')
                    f.write('y_axis'+'\t%.5f'*nbins%tuple(omega[a,b])+'\n')
                    f.close()
                rst_line = 'Dihedral CA %d CB %d CB %d CA %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,b+1,name,1.0,ASTEP)
                rst['omega'].append([a,b,p,rst_line])
        print("omega restraints: %d"%(len(rst['omega'])))
        ########################################################
        # theta: -pi..pi
        ########################################################
        prob = np.sum(theta[:,:,1:], axis=-1)
        i,j = np.where(prob>PCUT)
        prob = prob[i,j]
        theta = -np.log((theta+MEFF)/(theta[:,:,-1]+MEFF)[:,:,None])
        theta = np.concatenate([theta[:,:,-2:],theta[:,:,1:],theta[:,:,1:3]],axis=-1)
        for a,b,p in zip(i,j,prob):
            if b!=a:
                name=self.tmpdir.name+"/%d.%d_theta.txt"%(a+1,b+1)
                with open(name, "w") as f:
                    f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                    f.write('y_axis'+'\t%.3f'*nbins%tuple(theta[a,b])+'\n')
                    f.close()
                rst_line = 'Dihedral N %d CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,a+1,b+1,name,1.0,ASTEP)
                rst['theta'].append([a,b,p,rst_line])
        print("theta restraints: %d"%(len(rst['theta'])))
        ########################################################
        # phi: 0..pi
        ########################################################
        nbins = phi.shape[2]-1+4
        bins = np.linspace(-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
        prob = np.sum(phi[:,:,1:], axis=-1)
        i,j = np.where(prob>PCUT)
        prob = prob[i,j]
        phi = -np.log((phi+MEFF)/(phi[:,:,-1]+MEFF)[:,:,None])
        phi = np.concatenate([np.flip(phi[:,:,1:3],axis=-1),phi[:,:,1:],np.flip(phi[:,:,-2:],axis=-1)], axis=-1)
        for a,b,p in zip(i,j,prob):
            if b!=a:
                name=self.tmpdir.name+"/%d.%d_phi.txt"%(a+1,b+1)
                with open(name, "w") as f:
                    f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                    f.write('y_axis'+'\t%.3f'*nbins%tuple(phi[a,b])+'\n')
                    f.close()
                rst_line = 'Angle CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,name,1.0,ASTEP)
                rst['phi'].append([a,b,p,rst_line])
        print("phi restraints:   %d"%(len(rst['phi'])))
        return rst

    def select_rst(self, rst, dmin, dmax, intra=True, inter=True):
        array = []
        keys = list(rst.keys())
        pcut = self.par['PCUT']
        sep = self.par['LENSEQ1'] 
        seq = self.par['SEQ']
        for key in keys:
            for a,b,p,line in rst[key]:
                if abs(a-b)<dmin or abs(a-b)>=dmax: continue
                if key=='dist' and p<pcut: continue
                if key=='omega' and p<pcut+0.5: continue
                if key=='theta' and p<pcut+0.5: continue
                if key=='phi' and p<pcut+0.6: continue
                if not intra and np.sign(a-sep)==np.sign(b-sep): continue
                if not inter and np.sign(a-sep)!=np.sign(b-sep): continue
                if seq[a]=='G' or seq[b]=='G': continue
                array.append(line)
        return array

    def add_site(self, npz_file, array):
        npz = np.load(npz_file)
        dist = npz['dist'][:,:,5:]
        sep = self.par['LENSEQ1']
        pcut = self.par['PCUT']
        psumup = np.sum(dist, axis=-1)
        i,j = np.where(psumup>pcut+0.5)
        for a,b in zip(i,j):
            if a>b or abs(a-b)<=5 or np.sign(a-sep)==np.sign(b-sep): continue
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
        array = []
        for a in siteA:
            rstline = 'AmbiguousConstraint\n'
            for b in siteB:
                rstline += 'AtomPair CB {} CB {} FLAT_HARMONIC 10 1 4\n'.format(a, b)
            rstline += 'END\n'
            array.append(rstline)
        for a in siteB:
            rstline = 'AmbiguousConstraint\n'
            for b in siteA:
                rstline += 'AtomPair CB {} CB {} FLAT_HARMONIC 10 1 4\n'.format(a, b)
            rstline += 'END\n'
            array.append(rstline)
        array += ['AtomPair CB {} CB {} FLAT_HARMONIC 60 1 40'.format(a, b) for a in repA for b in repB]
        array += ['AtomPair CB {} CB {} FLAT_HARMONIC 60 1 50'.format(a, b) for a in repA for b in siteB]
        array += ['AtomPair CB {} CB {} FLAT_HARMONIC 60 1 50'.format(a, b) for a in repB for b in siteA]
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

        sf_fa = create_score_function('ref2015')
        sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)

        relax = rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sf_fa)
        relax.dualspace(True)
        relax.set_movemap(mmap)

        to_fullatom = SwitchResidueTypeSetMover('fa_standard')
        to_fullatom.apply(pose)

        siteA, repA = self.ISPRED_site(isfile1, 0)
        siteB, repB = self.ISPRED_site(isfile2, self.par['LENSEQ1'])
        array = self.format_ISPRED_rst(siteA, siteB, repA, repB)
        print ('Extracted {} constraints'.format(len(array)))
        self.apply_rst(pose, array, self.tmpdir.name+'/minimize.cst')

        print (sf_fa(pose))
        #emap = rosetta.core.scoring.EMapVector()
#        chain1 = range(pose.chain_begin(1), pose.chain_end(1))
#        chain2 = range(pose.chain_begin(2), pose.chain_end(2))
#        for r1 in chain1:
#            for r2 in chain2:
#                sf_fa.eval_ci_2b(pose.residue(r1), pose.residue(r2), pose, emap)
#                if emap[rosetta.core.scoring.site_constraint] != 0:
#                    print (emap[rosetta.core.scoring.site_constraint])

        relax.apply(pose)

#        for r1 in chain1: 
#            for r2 in chain2:
#                sf_fa.eval_ci_2b(pose.residue(r1), pose.residue(r2), pose, emap)
#        print (emap[rosetta.core.scoring.site_constraint])

        print (sf_fa(pose))
        pose.dump_pdb(ns.out+'_customr.pdb')

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



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

class constr_fold():

    def __init__(self, scriptdir, datadir):
        self.scriptdir = scriptdir.rstrip('/')+'/'
        self.datadir = datadir.rstrip('/')+'/'
        self.tmpdir = tempfile.TemporaryDirectory(prefix=datadir)

        with open(scriptdir+'/data/params.json') as jsonfile:
            self.par = json.load(jsonfile)

        init('-hb_cen_soft -relax:default_repeats 5 \
              -default_max_cycles 200 -out:level 100')

        sf = ScoreFunction()
        sf.add_weights_from_file(self.scriptdir + '/data/scorefxn.wts')
        sf1 = ScoreFunction()
        sf1.add_weights_from_file(self.scriptdir + '/data/scorefxn1.wts')
        sf_vdw = ScoreFunction()
        sf_vdw.add_weights_from_file(self.scriptdir + '/data/scorefxn_vdw.wts')
        sf_cart = ScoreFunction()
        sf_cart.add_weights_from_file(self.scriptdir + '/data/scorefxn_cart.wts')
        self.sfs = (sf, sf1, sf_vdw, sf_cart)


    def pose_from_fasta(self, fasta_file):
        self.sequence = ''
        for line in open(self.datadir+fasta_file): 
            if not line.startswith('>'): self.sequence += line.rstrip()
        pose = pose_from_sequence(self.sequence, 'centroid')
        return pose

    def join_fragment_poses(self, pose1, pose2, chain):
        pose1.delete_polymer_residue(len(pose1.sequence()))
        for i in range(1, len(pose2.sequence())):
            pose1.append_residue_by_bond(pose2.residue(i+1), True)
        pose1.pdb_info().set_chains(chain)
        renumber_pdbinfo_based_on_conf_chains(pose1)
        return pose1

    def mutate_gly(self, sequence, pose, mut):
        for i,a in enumerate(sequence):
            if a == 'G':
                mutator = rosetta.protocols.simple_moves.MutateResidue(i+1, mut.upper())
                mutator.apply(pose)

    def remove_clashes(self, pose, mover, rc_sf, thr=10, repeat=5):
        for _ in range(0, repeat):
            if float(rc_sf(pose)) < thr: break
            mover.apply(pose)

    def npz_selection(self, npz, chain):
        sep = self.par['LENSEQ1']
        assert chain in [1, 2], 'Chain must be 1 or 2!'
        if chain==1:
            npz = {'dist' : npz['dist'][0:sep, 0:sep],
                   'omega': npz['omega'][0:sep, 0:sep],
                   'theta': npz['theta'][0:sep, 0:sep],
                   'phi'  : npz['phi'][0:sep, 0:sep]}
        else:
            npz = {'dist' : npz['dist'][sep:, sep:],
                   'omega': npz['omega'][sep:, sep:],
                   'theta': npz['theta'][sep:, sep:],
                   'phi'  : npz['phi'][sep:, sep:]}
        return npz

    def format_rst(self, npz):
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

    def select_rst(self, rst, dmin, dmax, seq, intra=True, inter=True, limit=None):
        array = []
        keys = list(rst.keys())
        pcut = self.par['PCUT']
        sep = self.par['LENSEQ1'] 
        for key in keys:
            for a,b,p,line in rst[key]:
                if limit!=None: 
                    if a>=limit or b>=limit: continue
                if abs(a-b)<dmin or abs(a-b)>=dmax: continue
                if not intra and np.sign(a-sep)==np.sign(b-sep): continue
                if not inter and np.sign(a-sep)!=np.sign(b-sep): continue
                if key=='dist' and p<pcut+0.5: continue
                if key=='omega' and p<pcut+0.5: continue
                if key=='theta' and p<pcut+0.5: continue
                if key=='phi' and p<pcut+0.6: continue
                if seq[a]=='G' or seq[b]=='G': continue
                array.append(line)
        return array

    def apply_rst(self, pose, array, constraints_file):
        with open(constraints_file,'w') as f:
            for line in array: f.write(line+'\n')
        constraints = rosetta.protocols.constraint_movers.ConstraintSetMover()
        constraints.constraint_file(constraints_file)
        constraints.add_constraints(True)
        constraints.apply(pose)
        os.remove(constraints_file)

    def folding(self, pose, rst, mmap):

        min_mover = MinMover(mmap, self.sfs[0], 'lbfgs_armijo_nonmonotone', 0.0001, True)
        min_mover.max_iter(1000)
        min_mover1 = MinMover(mmap, self.sfs[1], 'lbfgs_armijo_nonmonotone', 0.0001, True)
        min_mover1.max_iter(1000)
        min_mover_vdw = MinMover(mmap, self.sfs[2], 'lbfgs_armijo_nonmonotone', 0.0001, True)
        min_mover_vdw.max_iter(500)
        min_mover_cart = MinMover(mmap, self.sfs[3], 'lbfgs_armijo_nonmonotone', 0.0001, True)
        min_mover_cart.max_iter(1000)
        min_mover_cart.cartesian(True)
        repeat_mover = RepeatMover(min_mover, 3)

        seqlen = len(pose.sequence())
        array = self.select_rst(rst, 1, seqlen, pose.sequence(), limit=seqlen)
        print (len(array))
        self.apply_rst(pose, array, self.tmpdir.name+'/minimize.cst')
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        self.remove_clashes(pose, min_mover1, self.sfs[2])

        return pose

    def incremental_folding(self, pose, rst, chain, init=100, step=40, overlap=60):
        fpose = None
        mmap = MoveMap()
        fragsidx = [[0,init]]
        fragsidx += [[f-(step+2),f] for f in range(init+step, len(pose.sequence()), step)]
        for frag in fragsidx: print (frag)
        frags = [pose.sequence()[s:e] for [s,e] in fragsidx]
        frags += [pose.sequence()[fragsidx[-1][-1]-2:]]
        if frags[-1] == '': frags = frags[:-1]
        #print (pose.sequence())
        #print (frags)

        for frag in frags:
            newfrag = pose_from_sequence(frag, 'centroid')
            if fpose == None: fpose = newfrag 
            else: fpose = self.join_fragment_poses(fpose, newfrag, chain)

            mmap.set_bb(False)
            end = len(fpose.sequence())
            #if frag == frags[-1]: end = len(fpose.sequence())
            #else: end = len(fpose.sequence())-overlap
            if frag == frags[0]: start = 1
            else: start = end-(step+overlap)
            mmap.set_bb_true_range(start,end)

            fpose.remove_constraints()
            fpose = self.folding(fpose, rst, mmap)
            fpose.dump_pdb('models/f{}{}_{}'.format(len(fpose.sequence()), chain, ns.out))
        return fpose

    def relax(self, pose, array, energy_func='ref2015', bb=True, chi=True, jump=True):
        sf_fa = create_score_function(energy_func)
        sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)
        sf_fa.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
        sf_fa.set_weight(rosetta.core.scoring.angle_constraint, 1)
    
        mmap = MoveMap()
        mmap.set_bb(bb)
        mmap.set_chi(chi)
        mmap.set_jump(jump)
    
        relax = rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sf_fa)
        relax.max_iter(200)
        relax.dualspace(True)
        relax.set_movemap(mmap)
    
        pose.remove_constraints()
        switch = SwitchResidueTypeSetMover("fa_standard")
        switch.apply(pose)
    
        self.par['PCUT'] = 0.15
        self.apply_rst(pose, array, self.tmpdir.name+'/minimize.cst')        
        relax.apply(pose)


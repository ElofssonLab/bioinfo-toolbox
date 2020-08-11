import sys,os,json
import tempfile
import numpy as np

from arguments import *
from utils_ros import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.protocols.rigid import *

os.environ["OPENBLAS_NUM_THREADS"] = "1"

def main():

    ########################################################
    # process inputs
    ########################################################

    # read params
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    #print (scriptdir)
    with open(scriptdir + '/data/params.json') as jsonfile:
        params = json.load(jsonfile)

    # get command line arguments
    args = get_args(params)
    #print(args)
    #print (params)
    params["interchain"]=args.interchain
    #print (params["interchain"])
    #sys.exit()
    # init PyRosetta
    init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')

    # Create temp folder to store all the restraints
    tmpdir = tempfile.TemporaryDirectory(prefix=args.wdir+'/')
    params['TDIR'] = tmpdir.name
    print('temp folder:     ', tmpdir.name)

    # read and process restraints & sequence
    npz = np.load(args.NPZ)
    seq1 = read_fasta(args.FASTA)
    seq2=read_fasta(args.FASTA2)
    seq=seq1+seq2
    L = len(seq)
    params['seq'] = seq
    params["seqlen1"]=len(seq1)
    params["seqlen2"]=len(seq2)
    params["seqlen"]=len(seq)
    rst = gen_rst(npz,tmpdir,params)
    #print (args.minprob)
    #sys.exit()
    #add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts) # Adding a weak flat harmonic to bring things together.
    #add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts) # Adding a weak flat harmonic to bring things together. 
    seq_polyala = 'A'*len(seq1+seq2) # Is this used ?

    #print (rst)
    #print (rst)
                       
    ########################################################
    # Scoring functions and movers
    ########################################################
    sf = ScoreFunction()
    sf.add_weights_from_file(scriptdir + '/data/scorefxn.wts')

    sf1 = ScoreFunction()
    sf1.add_weights_from_file(scriptdir + '/data/scorefxn1.wts')

    sf_vdw = ScoreFunction()
    sf_vdw.add_weights_from_file(scriptdir + '/data/scorefxn_vdw.wts')

    sf_cart = ScoreFunction()
    sf_cart.add_weights_from_file(scriptdir + '/data/scorefxn_cart.wts')

    mmap = MoveMap()
    mmap.set_bb(True)
    mmap.set_chi(False)
    mmap.set_jump(True)

    min_mover = MinMover(mmap, sf, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover.max_iter(1000)

    min_mover1 = MinMover(mmap, sf1, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover1.max_iter(1000)

    min_mover_vdw = MinMover(mmap, sf_vdw, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover_vdw.max_iter(500)

    min_mover_cart = MinMover(mmap, sf_cart, 'lbfgs_armijo_nonmonotone', 0.0001, True)
    min_mover_cart.max_iter(1000)
    min_mover_cart.cartesian(True)

    repeat_mover = RepeatMover(min_mover, 3)


    ########################################################
    # initialize pose
    ########################################################
    pose = pose_from_sequence(seq1, 'centroid' )
    pose2 = pose_from_sequence(seq2, 'centroid' )
    append_pose_to_pose(pose,pose2)
    # mutate GLY to ALA
    for i,a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'ALA')
            mutator.apply(pose)
            print('mutation: G%dA'%(i+1))

    set_random_dihedral(pose)
    remove_clash(sf_vdw, min_mover_vdw, pose)
    #pose.dump_pdb("starting.pdb")

    ########################################################  
    # minimization
    ########################################################

    
    #params["interchain"]=False
    
    if args.mode == 0:
        
        # short
        print('short')
        add_rst_chain2(pose, rst, 1, 12, params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

        # medium
        print('medium')
        add_rst_chain2(pose, rst, 12, 24, params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

        # long
        print('long')
        add_rst_chain2(pose, rst, 24, len(seq), params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

    elif args.mode == 1:

        # short + medium
        print('short + medium')
        add_rst_chain2(pose, rst, 3, 24, params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

        # long
        print('long')
        add_rst_chain2(pose, rst, 24, len(seq), params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

    elif args.mode == 2:

        # short + medium + long
        print('short + medium + long')
        add_rst_chain2(pose, rst, 1, len(seq), params)
        repeat_mover.apply(pose)
        min_mover_cart.apply(pose)
        remove_clash(sf_vdw, min_mover1, pose)

    #pose.dump_pdb("intraoptimized.pdb")

    #params["interchain"]=False
    #add_rst_chain2(pose, rst, 1, len(seq), params)
    #add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts)
    add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts) # Adding a weak flat harmonic to bring things together. 
    add_intra_rst(pose, rst, 1, len(seq), params)
    #repeat_mover.apply(pose)
    min_mover_cart.apply(pose)
    remove_clash(sf_vdw, min_mover1, pose)


    ## Now we add the interchaion contacts
    #    
    #mmap = MoveMap()
    #mmap.set_bb(False)
    #mmap.set_chi(False)
    #mmap.set_jump(False)
    #mmap.set_branches(False)
    #mmap.set_jump(1, True)
    #
    ## Mover from rigit body docking
    ##pert_mover = RigidBodyPerturbMover(jump_num, 3,8 )
    #
    #params["interchain"]=True
    #print('Moving the chains together')
    #add_rst_chain2(pose, rst, 1, len(seq), params)
    #repeat_mover.apply(pose)
    #min_mover_cart.apply(pose)
    #remove_clash(sf_vdw, min_mover1, pose)
    #
    #
    #minmover = MinMover()
    #minmover.movemap(movemap)
    #minmover.score_function(scorefxn) # use any scorefxn
    #minmover.apply(pose)
        
    # mutate ALA back to GLY
    for i,a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'GLY')
            mutator.apply(pose)
            print('mutation: A%dG'%(i+1))

    #pose.dump_pdb("interoptimized.pdb")
    ########################################################
    # full-atom refinement
    ########################################################

    if args.fastrelax == True:

        sf_fa = create_score_function('ref2015')
        sf_fa.set_weight(rosetta.core.scoring.atom_pair_constraint, 5)
        sf_fa.set_weight(rosetta.core.scoring.dihedral_constraint, 1)
        sf_fa.set_weight(rosetta.core.scoring.angle_constraint, 1)

        mmap = MoveMap()
        mmap.set_bb(True)
        mmap.set_chi(True)
        mmap.set_jump(True)

        relax = rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(sf_fa)
        relax.max_iter(200)
        relax.dualspace(True)
        relax.set_movemap(mmap)

        pose.remove_constraints()
        switch = SwitchResidueTypeSetMover("fa_standard")
        switch.apply(pose)

        print('relax...')
        params['PCUT'] = 0.15
        add_rst_chain2(pose, rst, 1, len(seq), params, True)
        relax.apply(pose)

    ########################################################
    # save final model
    ########################################################
    pose.dump_pdb(args.OUT)


if __name__ == '__main__':
    main()

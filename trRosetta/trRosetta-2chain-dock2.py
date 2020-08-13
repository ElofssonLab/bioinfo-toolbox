import sys,os,json
import tempfile
import numpy as np

from arguments import *
from utils_ros import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta.rosetta.protocols.docking import *
from pyrosetta.rosetta.protocols.simple_moves import *

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

    ##############################################################################
    # docking movers and energy functions
    ##############################################################################
    #scorefxn_docking_high = create_score_function('docking')
    #scorefxn_high_min = create_score_function('docking', 'docking_min')
    
    scorefxn_docking_low = create_score_function('interchain_cen')
    scorefxn_low_min = create_score_function('interchain_cen', 'docking_min')
    
    translation=3.0
    rotation = 8.0
    dock_jump=1
    dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
    spin = RigidBodySpinMover(dock_jump)
    slide_into_contact = DockingSlideIntoContact(dock_jump)
    
    movemap_docking = MoveMap()
    movemap_docking.set_jump(dock_jump, True)
    minmover_docking = MinMover()
    minmover_docking.movemap(movemap_docking)
    minmover_docking.score_function(scorefxn_low_min)
    #slide = DockingSlideIntoContact(dock_jump)


    ########################################################
    # initialize pose
    ########################################################
    pose = pose_from_sequence(seq1, 'centroid' )
    pose2 = pose_from_sequence(seq2, 'centroid' )
    append_pose_to_pose(pose,pose2)
    partners="A_B"
    dock_jump=1
    # mutate GLY to ALA
    for i,a in enumerate(seq):
        if a == 'G':
            mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'ALA')
            mutator.apply(pose)
            print('mutation: G%dA'%(i+1))
            
    # The order of these two commands matters, to make the proteins separated or not        
    set_random_dihedral(pose)
    setup_foldtree(pose, partners, Vector1([dock_jump]))

    # Just to ensure that the two proteins starts separated.
    mindist=100
    #print ("Rotation 1: ",pose.jump(dock_jump).get_rotation())
    print ("Translation 1: ",pose.jump(dock_jump).get_translation())
    if sum(pose.jump(dock_jump).get_translation()) < mindist:
        initial_mover = RigidBodyPerturbMover(dock_jump, 0, 2*mindist)
        #print ("Rotation 1b: ",pose.jump(dock_jump).get_rotation())
        print ("Translation 1b: ",pose.jump(dock_jump).get_translation())
        
    remove_clash(sf_vdw, min_mover_vdw, pose)
    if args.saveintermediate:
        pose.dump_pdb(args.OUT+"-starting.pdb")
    starting_pose=Pose()
    starting_pose.assign(pose)
    ########################################################  
    # minimization
    ########################################################

    # Alternative methods (that might work) (and flags)
    # Initital step
    # 1. Only predicted interactions using fade (will results in proteins separated when no interactions present)
    # 2. Add weak harmonic force to bring proteins together on all
    # 2a. Interacting pairs Fade function only
    # 2b. Interacting pairs with harmonic
    # 2c. Interacting pairs with harmonic + all residues (weak harmonic)
    # 2d. Predicted interaction surfaces  (with stronger potential?)
    # We also have options "allcontacts" for 2a and 2d and "interchain" for all versions
    # 3. Docking methods (i.e. here we add interaction potentials)
    # 3a. No docking
    # 3b. Old method (i.e. pyrosetta protocol a second run)
    # 3c. Minmover docking
    # 3d. dock_prot docking
    
    
    #params["interchain"]=False
    
    # Flags (lets use new that overrides some old ones)
    #
    #if args.initialmode == A: # No harmonic intrachanin attraction
    #    params["interchain"]=True
    if args.initialmode == "B": # (Default)
        #params["interchain"]=True
        add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=False)
    elif args.initialmode == "C":
        #params["interchain"]=True
        add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=True)
    elif args.initialmode == "D": 
        #params["interchain"]=True
        #add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=True)
        add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts)

    
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

    #print ("Rotation 2: ",pose.jump(dock_jump).get_rotation())
    print ("Translation 2: ",pose.jump(dock_jump).get_translation())

    if args.saveintermediate:
        pose.dump_pdb(args.OUT+"-step1.pdb")

    # Now we do dosking

    # Provide a random orientation of the two chains
    randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump,    partner_upstream)
    randomize_downstream = RigidBodyRandomizeMover(pose, dock_jump,   partner_downstream)

    random_pose=Pose()
    random_pose.assign(pose)
    
    
    mmap = MoveMap()
    mmap.set_bb(False)
    mmap.set_chi(False)
    mmap.set_jump(True)
    
    
    perturb = SequenceMover()
    perturb.add_mover(randomize_upstream)
    perturb.add_mover(randomize_downstream)
    perturb.add_mover(dock_pert)
    perturb.add_mover(spin)
    perturb.add_mover(slide_into_contact)
    perturb.add_mover(minmover_docking)
    
    
    dock_prot = DockingProtocol()    # contains many docking functions
    dock_prot.set_movable_jumps(Vector1([1]))    # set the jump to jump 1
    dock_prot.set_lowres_scorefxn(sf_vdw) # This should probably be sf1 or something like that
    dock_prot.set_low_res_protocol_only # This should probably be sf1 or something like that
    #dock_prot.set_highres_scorefxn(scorefxn_high_min)
    
    
    
    
    #repeat_dock = RepeatMover(min_mover1, 3)
    minmover_docking = MinMover()
    minmover_docking.movemap(movemap_docking)
    minmover_docking.score_function(sf_vdw)
    
    minmover_docking_sf1 = MinMover()
    minmover_docking_sf1.movemap(movemap_docking)
    minmover_docking_sf1.score_function(sf1)
    
    repeat_dock = SequenceMover()
    #repeat_dock.add_mover(randomize_upstream)
    #repeat_dock.add_mover(randomize_downstream)
    repeat_dock.add_mover(dock_pert)
    repeat_dock.add_mover(slide_into_contact)
    repeat_dock.add_mover(spin)
    repeat_dock.add_mover(slide_into_contact)
    repeat_dock.add_mover(spin)
    repeat_dock.add_mover(minmover_docking)
    #repeat_dock.add_mover(minmover_docking_sf1)
    
    
    step1_pose=Pose()
    step1_pose.assign(pose)
    


    params["interchain"]=True  # Always true here

    #add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts)
    #\ Adding a weak flat harmonic to bring things together.


    #params["interchain"]=False
    #add_rst_chain2(pose, rst, 1, len(seq), params)
    #add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts)
    #add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts) # Adding a weak flat harmonic to bring things together. 
    #add_intra_rst(pose, rst, 1, len(seq), params)
    #repeat_mover.apply(pose)


    # Sometimes we need to add the intra restraints again
    if args.initialmode == "A": # Should we also add it for "D"
        add_intra_rst(pose, rst, 1, len(seq), params)

    #  Dockingmode "A" is no docking  (# in principle we could try all modes at the same time
    if "B" in args.dockingmode: # default protocol
        dock_pose=Pose()
        dock_pose.assign(pose)
        min_mover_cart.apply(dock_pose)
        if args.saveintermediate:
            dock_pose.dump_pdb(args.OUT+"-dockB1.pdb")
        remove_clash(sf_vdw, min_mover1, dock_pose)
        if args.saveintermediate:
            dock_pose.dump_pdb(args.OUT+"-dockB2.pdb")
            #pose.dump_pdb("intraoptimized.pdb")
        print ("Rotation 3B: ",pose.jump(dock_jump).get_rotation())
        print ("Translation 3B: ",pose.jump(dock_jump).get_translation())
    if "C" in args.dockingmode: # Docking protocol using minmover
        dock_pose=Pose()
        dock_pose.assign(pose)
        repeat_dock.apply(dock_pose)
        if args.saveintermediate:
            dock_pose.dump_pdb(args.OUT+"-dockC1.pdb")
        remove_clash(sf_vdw, min_mover1, dock_pose)
        if args.saveintermediate:
            dock_pose.dump_pdb(args.OUT+"-dockC2.pdb")
        print ("Rotation 3C: ",pose.jump(dock_jump).get_rotation())
        print ("Translation 3C: ",pose.jump(dock_jump).get_translation())
    if "D" in args.dockingmode: # Docking protocol using minmover
        dock_pose=Pose()
        dock_pose.assign(pose)
        dock_prot.apply(dock_pose)
        if args.saveintermediate:
            dock_pose.dump_pdb(args.OUT+"-dockD1.pdb")
        remove_clash(sf_vdw, min_mover1, pose)
        if args.saveintermediate:
            dock_pose.dump_pdb(args.OUT+"-dockD2.pdb")
        print ("Rotation 3D: ",pose.jump(dock_jump).get_rotation())
        print ("Translation 3D: ",pose.jump(dock_jump).get_translation())


        
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

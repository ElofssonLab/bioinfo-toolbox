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

model="1ay7"

########################################################
# process inputs
########################################################
def get_args_default(params):
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-NPZ", type=str,default="pydocktest/"+model+"_u1_A-"+model+"_u2_A.npz", help="input distograms and anglegrams (NN predictions)")
    parser.add_argument("-FASTA", type=str,default="pydocktest/"+model+"_u1_A.fasta", help="input sequence")
    parser.add_argument("-FASTA2", type=str,default="pydocktest/"+model+"_u2_A.fasta", help="input sequence 2")
    parser.add_argument("-OUT", type=str,default="pydocktest/"+model+"_u1_A-"+model+"_u2_A.pdb", help="output model (in PDB format)")

    
    parser.add_argument('-minprob', type=float, dest='minprob', default=0.5, help='min probability of distance restraints for inter-chain flat harminic')
    parser.add_argument('-intradist', type=float, dest='intradist', default=15., help='The distance for the flat part in  restraints for inter-chain flat harmonic')
    parser.add_argument('-intrasd', type=float, dest='intrasd', default=20., help='SD for  restraints for inter-chain flat harminic')
    parser.add_argument('-allintra',  dest='allcontacts', help='Use in additiona also a weak interacting potential for all intrachain contacts',action='store_true')
    parser.add_argument('-pd', type=float, dest='pcut', default=params['PCUT'], help='min probability of distance restraints')
    parser.add_argument('-m', type=int, dest='mode', default=2, choices=[0,1,2], help='0: sh+m+l, 1: (sh+m)+l, 2: (sh+m+l)')
    parser.add_argument('-w', type=str, dest='wdir', default=params['WDIR'], help='folder to store temp files')
    parser.add_argument('-n', type=int, dest='steps', default=1000, help='number of minimization steps')
    parser.add_argument('--orient', dest='use_orient', action='store_true', help='use orientations')
    parser.add_argument('--nointerchain', dest='interchain', action='store_false', help='',default=True)
    parser.add_argument('--no-orient', dest='use_orient', action='store_false')
    parser.add_argument('--fastrelax', dest='fastrelax', action='store_true', help='perform FastRelax')
    parser.add_argument('--no-fastrelax', dest='fastrelax', action='store_false')
    parser.set_defaults(use_orient=True)
    parser.set_defaults(fastrelax=True)

    args = parser.parse_args()

    params['PCUT'] = args.pcut
    params['USE_ORIENT'] = args.use_orient

    return args

# read params
scriptdir = "/home/arnee/git/bioinfo-toolbox/trRosetta/"
#print (scriptdir)
with open(scriptdir + '/data/params.json') as jsonfile:
    params = json.load(jsonfile)

# get command line arguments
args = get_args_default(params)
#print(args)

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
#add_intrachain_rst(npz,rst,tmpdir,params,minprob=0) # Adding a weak flat harmonic to bring things together. 
#seq_polyala = 'A'*len(seq1+seq2) # Is this used ?

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
dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
spin = RigidBodySpinMover(dock_jump)
dock_jump=1
slide_into_contact = DockingSlideIntoContact(dock_jump)

movemap_docking = MoveMap()
movemap_docking.set_jump(dock_jump, True)
minmover_docking = MinMover()
minmover_docking.movemap(movemap_docking)
minmover_docking.score_function(scorefxn_low_min)



slide = DockingSlideIntoContact(1)

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

    
remove_clash(sf_vdw, min_mover_vdw, pose)
pose.dump_pdb("starting.pdb")

########################################################  
# minimization
########################################################


# Alternative methods (that might work)
# 1. use interchina=False first, then docking with interachain=True 
# 2. use interchina=False first, then docking with interachain=True + add_instrachain_rst
# 3. use interchina=False first, then docking with interachain=True + add_instrachain_rst + allintra
# 4. use interchina=True (works
# 5. use interchina=True + add_instrachain_rst + allintra - no docking
# 6. use interchina=True + add_instrachain_rst + allintra + docking


params["interchain"]=False
add_intrachain_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts)
# Adding a weak flat harmonic to bring things together. 
starting_pose=Pose()
starting_pose.assign(pose)
# short + medium + long
print('short + medium + long')
add_rst_chain2(pose, rst, 1, len(seq), params)
repeat_mover.apply(pose)
min_mover_cart.apply(pose)
remove_clash(sf_vdw, min_mover1, pose)

pose.dump_pdb("step1.pdb")

step1_pose=Pose()
step1_pose.assign(pose)


# Provide a random orientation of the two chains
randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump,    partner_upstream)
randomize_downstream = RigidBodyRandomizeMover(pose, dock_jump,   partner_downstream)

#randomorient = SequenceMover()
#randomorient.add_mover(randomize_upstream)
#randomorient.add_mover(dock_pert)
#randomorient.apply(pose)
#remove_clash(sf_vdw, min_mover1, pose)
#pose.dump_pdb("step2.pdb")
#random_pose=Pose()
#random_pose.assign(pose)


params["interchain"]=True
add_interacton_areas_rst(npz,rst,tmpdir,params,minprob=args.minprob,UB=args.intradist,D=args.intrasd,allcontacts=args.allcontacts)
#\ Adding a weak flat harmonic to bring things together.
add_intra_rst(pose, rst, 1, len(seq), params)
#repeat_mover.apply(pose)
#min_mover_cart.apply(pose)
#remove_clash(sf_vdw, min_mover1, pose)
#pose.dump_pdb("step3.pdb")

# We use energey function sf_vdw (or sf1)

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
repeat_dock1 = SequenceMover()
repeat_dock2 = SequenceMover()
repeat_dock3 = SequenceMover()
repeat_dock4 = SequenceMover()
repeat_dock5 = SequenceMover()
repeat_dock6 = SequenceMover()
repeat_dock7 = SequenceMover()
#repeat_dock.add_mover(randomize_upstream)
#repeat_dock.add_mover(randomize_downstream)
repeat_dock.add_mover(dock_pert)
repeat_dock.add_mover(slide_into_contact)
repeat_dock.add_mover(spin)
repeat_dock.add_mover(slide_into_contact)
repeat_dock.add_mover(spin)
repeat_dock.add_mover(minmover_docking)
repeat_dock1.add_mover(dock_pert)
repeat_dock2.add_mover(slide_into_contact)
repeat_dock3.add_mover(spin)
repeat_dock4.add_mover(slide_into_contact)
repeat_dock5.add_mover(spin)
repeat_dock6.add_mover(minmover_docking)
repeat_dock7.add_mover(minmover_docking_sf1)


pose=Pose()
pose.assign(random_pose)
    
#repeat_dock.apply(pose)   # Brings them closer but random orientation
repeat_dock1.apply(pose)
repeat_dock2.apply(pose)
repeat_dock3.apply(pose)
repeat_dock4.apply(pose)
repeat_dock5.apply(pose)
repeat_dock6.apply(pose)
repeat_dock7.apply(pose)
repeat_dock.apply(pose)
repeat_dock.apply(pose)
repeat_dock.apply(pose)
pose.dump_pdb("minmover.pdb")
#min_mover.apply(pose)   # This unfolds the proteisn
#pose.dump_pdb("minmover-docking2.pdb")
remove_clash(sf_vdw, min_mover1, pose)
pose.dump_pdb("minmover2.pdb")  # This actually moves the packing and can fix the packing.. KEEP
minmover_pose=Pose()
minmover_pose.assign(pose)


# Extra dockng (test
pose=Pose()
pose.assign(random_pose)

dock_prot.apply(pose)
pose.dump_pdb("minmover-dock1.pdb") # For some reason this is full atom
#min_mover1.apply(pose)
#remove_clash(sf_vdw, min_mover1, pose)
#pose.dump_pdb("minmover-dock2.pdb") 


#jobs=1
#jd = PyJobDistributor(job_output, jobs, scorefxn_low_min)
#while not jd.job_complete:
#        # a. set necessary variables for this trajectory
#        # -reset the test pose to original (centroid) structure
#        test_pose.assign(pose)
#        # -change the pose name, for pretty output to PyMOL
#        counter += 1
#        test_pose.pdb_info().name(job_output + '_' + str(counter))
#
#        # b. perturb the structure for this trajectory
#        perturb.apply(test_pose)
#
#        # c. perform docking
#        dock_prot.apply(test_pose)
#        #### alternate application of the DockingProtocol pieces
#        #docking_low.apply(test_pose)
#        #docking_high.apply(test_pose)
#
#        # d. output the decoy structure
#        #to_fullatom.apply(test_pose)    # ensure the output is fullatom
#        # to PyMOL
#        test_pose.pdb_info().name(job_output + '_' + str( counter ) + '_fa')
#        # to a PDB file
#        jd.output_decoy(test_pose)
    

# Docking

# mutate ALA back to GLY
for i,a in enumerate(seq):
    if a == 'G':
        mutator = rosetta.protocols.simple_moves.MutateResidue(i+1,'GLY')
        mutator.apply(pose)
        print('mutation: A%dG'%(i+1))

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



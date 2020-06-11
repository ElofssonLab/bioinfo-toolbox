import sys,os,json
sys.path.append("/home/arnee/git/bioinfo-toolbox/trRosetta/")
import tempfile
import numpy as np

from arguments import *
from utils_ros import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.protocols.rigid import *

# read params
scriptdir = os.path.dirname("/home/arnee/git/bioinfo-toolbox/trRosetta/")
with open(scriptdir + '/data/params.json') as jsonfile:
    params = json.load(jsonfile)

# get command line arguments
#args = get_args(params)

# init PyRosetta
init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')


########################################################
# Scoring functions and movers
########################################################
# Create temp folder to store all the restraints
tmpdir = tempfile.TemporaryDirectory(prefix='/tmp/')
params['TDIR'] = tmpdir.name
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

#------------------------------------------------------------------------

i="nsp8"
j="nsp7"
# read and process restraints & sequence
npz = np.load(i+"-"+j+"-nosep.npz")
seq1 = read_fasta(i+".fa")
seq2=read_fasta(j+".fa")
seq=seq1+seq2
L = len(seq)

params['seq'] = seq
params["seqlen1"]=len(seq1)
params["seqlen2"]=len(seq2)
params["seqlen"]=len(seq)
rst = gen_rst(npz,tmpdir,params)
#seq_polyala = 'A'*len(seq1+seq2) # Is this used ?
add_intrachain_rst(rst,tmpdir,params)


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
pose.dump_pdb("starting.pdb")
params["USE_ORIENT"]=True
params["interchain"]=True
print(pose)
mmap = MoveMap()                            
mmap.set_bb(True)                                          
mmap.set_chi(True)      
mmap.set_jump(True) 
print (mmap)

#
#pose=pose_from_pdb("starting.pdb")

########################################################
# minimization
########################################################


def extra_Stuff():
    params["interchain"]=False
    params["USE_ORIENT"]=True    
    # short + medium + long
    print('short + medium + long')
    add_rst_chain2(pose, rst, 1, len(seq), params)
    repeat_mover.apply(pose)
    min_mover_cart.apply(pose)
    remove_clash(sf_vdw, min_mover1, pose)
    
    
    pose.dump_pdb("intraoptimized.pdb")
    
    # Now we add the interchaion contacts
    
    mmap = MoveMap()
    mmap.set_bb(False)
    mmap.set_chi(False)
    mmap.set_jump(False)
    mmap.set_branches(False)
    mmap.set_jump(1, True)
    
    print (pose)
    print (pose.jump(1))
    
    # Mover from rigit body docking
    
    print('Moving the chains together')
    params["interchain"]=True
    add_rst_chain2(pose, rst, 1, len(seq), params)
    pert_mover = RigidBodyPerturbMover(1, 3,8 )
    #repeat_mover.apply(pose)
    min_mover_cart.apply(pose)
    remove_clash(sf_vdw, min_mover1, pose)
    
    
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

    pose.dump_pdb("interoptimized.pdb")

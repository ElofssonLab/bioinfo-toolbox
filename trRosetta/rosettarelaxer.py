import sys,os,json
import tempfile
import numpy as np

from arguments import *
from utils_ros import *
from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.protocols.rigid import *

init('-hb_cen_soft -relax:default_repeats 5 -default_max_cycles 200 -out:level 100')


pose = pose_from_pdb(sys.argv[1])

scriptdir = os.path.dirname(os.path.realpath(__file__))
with open('/home/gabriele/Desktop/programs/trRosetta/data/params.json') as jsonfile:
    params = json.load(jsonfile)

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
relax.apply(pose)

pose.dump_pdb(sys.argv[1].rstrip('.pdb')+'_relaxed.pdb')

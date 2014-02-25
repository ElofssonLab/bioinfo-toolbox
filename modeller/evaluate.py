import sys

from modeller import *
from modeller.scripts import complete_pdb

model = sys.argv[1]
#scoref = sys.argv[2]

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
#mdl = complete_pdb(env, 'TvLDH.B99990002.pdb')
mdl = complete_pdb(env, model)

# Assess with DOPE:
s = selection(mdl)   # all atom selection
#s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',
#              normalize_profile=True, smoothing_window=15)
#s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=scoref,
#              normalize_profile=True, smoothing_window=15)
s.assess_dope()

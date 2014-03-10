import sys

from modeller import *
from modeller.scripts import complete_pdb

model_lst = []
score_lst = []

if len(sys.argv) == 2:
    model_lst.append(sys.argv[1])
elif len(sys.argv) > 2:
    model_lst += sys.argv[1:]
#scoref = sys.argv[2]

for model in model_lst:
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
    dope = s.assess_dope()
    score_lst.append((dope, model))

score_lst.sort(key=lambda itm: itm[0])

for sc in score_lst:
    dope = sc[0]
    model = sc[1]
    print "DOPE_SCORE: " + str(dope) + " " + model

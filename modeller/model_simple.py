import sys

from modeller import *
from modeller.automodel import *

alnf = sys.argv[1]
id = sys.argv[2]
templ_id = sys.argv[3]

env = environ()
a = automodel(env, alnfile=alnf,
              knowns=templ_id, sequence=id,
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 100
a.make()

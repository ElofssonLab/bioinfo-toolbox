import sys

from modeller import *

query_id = sys.argv[1]
query_seq = sys.argv[2]
templ_id = sys.argv[3]
templ_pdb = sys.argv[4]
ali_file = sys.argv[5]

env = environ()
aln = alignment(env)
#mdl = model(env, file=templ_id, model_segment=('FIRST:A','LAST:A'))
mdl = model(env, file=templ_id, model_segment=('FIRST:D','LAST:D'))
aln.append_model(mdl, align_codes=templ_id, atom_files=templ_pdb)
aln.append(file=query_seq, align_codes=query_id)
aln.align2d()
aln.write(file=ali_file, alignment_format='PIR')

import sys

from modeller import *
from modeller.automodel import *
from modeller.optimizers import conjugate_gradients, molecular_dynamics, quasi_newton, actions

sys.path.append('/home/x_mirmi/bioinfo-toolbox')
from parsing import parse_contacts
from parsing import parse_psipred

query_id = sys.argv[1]
query_seq = sys.argv[2]
contact_filename = sys.argv[3]
psipred_filename = sys.argv[4]
factor = float(sys.argv[5])

env = environ()
env.io.atom_files_directory = '../atom_files'
env.edat.dynamic_sphere = True
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')


class MyFade(forms.restraint_form):
    """An implementation of Rosetta's FADE function"""

    rt = 0.5900991    # RT at 297.15K, in kcal/mol

    def __init__(self, group, feature, cutoff_lower, cutoff_upper, fade_zone, well_depth):
        forms.restraint_form.__init__(self, group, feature, 0, (cutoff_lower, cutoff_upper, fade_zone, well_depth))

    def eval(self, feats, iftyp, modal, param, deriv):
        (cutoff_lower, cutoff_upper, fade_zone, well_depth) = param
        z = feats[0]
        val = 1.0
        fderv = 0.0
        if z < cutoff_lower or z > cutoff_upper:
            val = 0.0
            fderv = 0.0
        elif z < (cutoff_lower + fade_zone):
            b = -1.0 * (z - (cutoff_lower + fade_zone)) / fade_zone
            b2 = b * b
            b3 = b2 * b
            val = 2 * b3 - 3 * b2 + 1
            fderv = -1.0 * (6 * b2 - 6 * b) / fade_zone
        elif z > (cutoff_upper - fade_zone):
            b = (z - (cutoff_upper - fade_zone)) / fade_zone
            b2 = b * b
            b3 = b2 * b
            val = 2 * b3 - 3 * b2 + 1
            fderv = (6 * b2 - 6 * b) / fade_zone
        val *= well_depth
        fderv *= well_depth
        if deriv:
            return val, [fderv]
        else:
            return val

    def vmin(self, feats, iftyp, modal, param):
        viol = 0.0
        (cutoff_lower, cutoff_upper, fade_zone, well_depth) = param
        z = feats[0]
        if z < (cutoff_lower + fade_zone):
            viol = abs(z - (cutoff_lower + fade_zone))
        if z > (cutoff_upper - fade_zone):
            viol = abs(z - (cutoff_upper - fade_zone))
        return viol

    def rvmin(self, feats, iftyp, modal, param):
        viol = 0.0
        (cutoff_lower, cutoff_upper, fade_zone, well_depth) = param
        z = feats[0]
        if z < (cutoff_lower + fade_zone):
            viol = abs(z - (cutoff_lower + fade_zone))
        if z > (cutoff_upper - fade_zone):
            viol = abs(z - (cutoff_upper - fade_zone))
        return viol

    def min_mean(self, feats, iftyp, modal, param):
        (cutoff_lower, cutoff_upper, fade_zone, well_depth) = param
        return [(cutoff_upper - abs(cutoff_lower)) / 2]

    def get_range(self, iftyp, modal, param, spline_range):
        (cutoff_lower, cutoff_upper, fade_zone, well_depth) = param
        #return (cutoff_lower, cutoff_upper)
        return (float('-inf'), float('inf'))

    # There is only one minimum, so the 'heavy' mean is the same as the 'min'
    vheavy = vmin
    rvheavy = rvmin
    heavy_mean = min_mean



# Create a new empty alignment and model:
aln = alignment(env)
mdl = model(env)

# Build a model from one-letter query_ids, and write to a PDB file:
aln.append(file=query_seq, align_codes=query_id)
mdl.generate_topology(aln[query_id])
mdl.build(initialize_xyz=True, build_method='INTERNAL_COORDINATES')
#mdl.build(initialize_xyz=True, build_method='3D_INTERPOLATION')
#mdl.build_sequence(aln)

# Add PDB remarks for readability
mdl.remark = """REMARK   4 Extended-chain model of 1fas
REMARK   4 Built from internal coordinates only"""

# Select all atoms:
atmsel = selection(mdl)

# Generate the restraints:
rsr = mdl.restraints
rsr.clear()
rsr.make(atmsel, restraint_type='stereo', spline_on_site=False)
ss_seq = parse_psipred.horizontal(open(psipred_filename, 'r'))
i = 1
"""
for ss in ss_seq:
    print ss
    if ss == 'H':
        rsr.add(secondary_structure.alpha(mdl.residue_range('%d:' % i, '%d:' % i)))
        #rsr.make(aln, restraint_type='ALPHA', residue_ids=('%d' % i, '%d' % i), spline_on_site=False)
    elif ss == 'E':
        rsr.add(secondary_structure.strand(mdl.residue_range('%d:' % i, '%d:' % i)))
    i += 1
"""
atmsel = selection(mdl)
mpdf = atmsel.energy()

# Create optimizer objects and set defaults for all further optimizations
cg = conjugate_gradients(output='REPORT')
md = molecular_dynamics(output='REPORT')

# Open a file to get basic stats on each optimization
trcfil = file(query_id+'.D00000001', 'w')

# Run CG on the all-atom selection; write stats every 5 steps
cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil), min_atom_shift=0.01)
# Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
# 10 steps during the run, and write stats every 10 steps
md.optimize(atmsel, temperature=300, max_iterations=50,
            actions=[actions.write_structure(10, query_id+'.D9998%04d.pdb'),
                     actions.trace(10, trcfil)])
# Finish off with some more CG, and write stats every 5 steps
cg.optimize(atmsel, max_iterations=20,
            actions=[actions.trace(5, trcfil)])

mpdf = atmsel.energy()

mdl.write(file=query_id+'.D00000001.pdb')

contacts = parse_contacts.parse(open(contact_filename, 'r'))
count = 0
seq_len = len(aln[query_id])
for (score, i, j) in contacts:
    rsr.add(forms.gaussian(group=physical.xy_distance,
                feature=features.distance(mdl.atoms['CA:%d' % i],
                                          mdl.atoms['CA:%d' % j]),
                mean=10.0, stdev=2))
    #rsr.add(MyFade(group=physical.xy_distance,
    #            feature=features.distance(mdl.atoms['CA:%d' % i],
    #                                      mdl.atoms['CA:%d' % j]),
    #            cutoff_lower=-100, cutoff_upper=100, fade_zone=92, well_depth=-150))
    is_gly_a = aln[query_id].residues[i-1].code == 'G'
    is_gly_b = aln[query_id].residues[j-1].code == 'G'
    """
    if is_gly_a and (not is_gly_b):
        rsr.add(forms.gaussian(group=physical.xy_distance,
                    feature=features.distance(mdl.atoms['CA:%d' % i],
                                              mdl.atoms['CB:%d' % j]),
                    mean=5.0, stdev=0.1))
    elif (not is_gly_a) and is_gly_b:
        rsr.add(forms.gaussian(group=physical.xy_distance,
                    feature=features.distance(mdl.atoms['CB:%d' % i],
                                              mdl.atoms['CA:%d' % j]),
                    mean=5.0, stdev=0.1))
    elif  (not is_gly_a) and (not is_gly_b):
        rsr.add(forms.gaussian(group=physical.xy_distance,
                    feature=features.distance(mdl.atoms['CB:%d' % i],
                                              mdl.atoms['CB:%d' % j]),
                    mean=5.0, stdev=0.1))
    """
    #rsr.add(forms.gaussian(group=physical.xy_distance,
    #            feature=features.distance(mdl.residues['%d' % i],
    #                                      mdl.residues['%d' % j]),
    #            mean=10.0, stdev=0.1))
    if count > seq_len * factor:
        break 
    count += 1

rsr.write(file=query_id+'.rsr')

atmsel = selection(mdl)
mpdf = atmsel.energy()

# Create optimizer objects and set defaults for all further optimizations
cg = conjugate_gradients(output='REPORT')
qn = quasi_newton(output='REPORT')
md = molecular_dynamics(output='REPORT')

# Open a file to get basic stats on each optimization
trcfil = file(query_id+'.D00000002', 'w')

# Run CG on the all-atom selection; write stats every 5 steps
#cg.optimize(atmsel, max_iterations=200, actions=actions.trace(5, trcfil), min_atom_shift=0.01)
qn.optimize(atmsel, max_iterations=200, actions=actions.trace(5, trcfil), min_atom_shift=0.01)
# Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
# 10 steps during the run, and write stats every 10 steps
md.optimize(atmsel, temperature=300, max_iterations=50,
            actions=[actions.write_structure(10, query_id+'.D9999%04d.pdb'),
                     actions.trace(10, trcfil)])
# Finish off with some more CG, and write stats every 5 steps
#cg.optimize(atmsel, max_iterations=20, actions=[actions.trace(5, trcfil)], min_atom_shift=0.01)
qn.optimize(atmsel, max_iterations=20, actions=[actions.trace(5, trcfil)], min_atom_shift=0.01)

mpdf = atmsel.energy()

mdl.write(file=query_id+'.D00000002.pdb')
"""

a = automodel(env,
              alnfile  = query_seq,     # alignment filename
              knowns   = query_id,              # codes of the templates
              sequence = query_id,              # code of the target
              csrfile  = query_id+'.rsr')            # use 'my' restraints file
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do homology modeling

"""
print count
print seq_len
print factor


import sys, os

from os.path import expanduser 
home = expanduser("~")         
sys.path.append(home + '/bioinfo-toolbox') 
 
from parsing import parse_contacts          
from parsing import parse_psipred           
from parsing import parse_fasta  


def is_beta(res_i, ss):
    result = False
    for offset in range(-1, 2):
    #for offset in range(0):
        try:
            result = result or ss[res_i + offset] == 'E'
        except IndexError:
            # sequence out of bounds: okay, pass it
            pass
    return result

# command line input
seqfile_name = sys.argv[1]
infile_name = sys.argv[2]
min_dist = int(sys.argv[3])

# set length of a single repeat unit:
#rep_len = 33

# use L * factor highest scoring constraints
factor = float(sys.argv[4])

if len(sys.argv) == 6:
    psipred_filename = sys.argv[5]
else:
    psipred_filename = ''

# scale value x from [min_x, max_x] to [0,1]
# NOT USED
def scale(x, min_x, max_x):
    range_x = (max_x - min_x)
    #print range_x
    return (x - min_x) / range_x


# read sequence needed for checking if res = glycine
# to assign CA/CB
seq = ''
seqfile = open(seqfile_name, 'r')
for line in seqfile:
    if line[0] != '>':
        seq += line.strip()
seqfile.close()
seq_len = len(seq)


# guessing separator of constraint file
test_line = open(infile_name,'r').readline()
if len(test_line.split(',')) != 1:
    sep = ','
elif len(test_line.split(' ')) != 1:
    sep = ' '
else:
    sep = '\t'


# read constraint file
infile = open(infile_name, 'r')
old_constraints = []
for line in infile:
    if not line.split():
        continue
    old_constraints.append(tuple(line.strip().split(sep)))
    
if psipred_filename:
    ss = parse_psipred.horizontal(open(psipred_filename, 'r'))
else:
    ss = ''


# reformat to rosetta constraints
old_constraints.sort(key=lambda x: float(x[-1]), reverse=True)
min_score = float(old_constraints[-1][-1])
max_score = float(old_constraints[0][-1])
rosetta_lines = []
count = 0
for constr in old_constraints:
    res1 = int(constr[0])
    res2 = int(constr[1])
    if abs(res2 - res1) >= min_dist: # and abs(res2 - res1) < rep_len * 1.5:
        atm1 = 'CB'
        if seq[res1 - 1] == 'G':
            atm1 = 'CA'
        atm2 = 'CB'
        if seq[res2 - 1] == 'G':
            atm2 = 'CA'
        #score = scale(float(constr[2]), min_score, max_score)
        #score = float(constr[-1])
        #score = (score * -20.0)# - (1/factor)
        score = -15.0
        if is_beta(res1 - 1, ss) and is_beta(res2 - 1, ss):
            print "Add beta-beta constraint at %s,%s" % (res1, res2)
            rosetta_lines.append('AtomPair %s %d %s %d FADE -10 19 10 %.2f 0' % (atm1, res1, atm2, res2, round(score, 2)*3))
        else:
            rosetta_lines.append('AtomPair %s %d %s %d FADE -10 19 10 %.2f 0' % (atm1, res1, atm2, res2, round(score, 2)))
        #rosetta_lines.append('AtomPair CA %d CA %d FADE -10 19 10 %.2f 0' % (res1, res2, round(score, 2)))
        #rosetta_lines.append('AtomPair %s %d %s %d BOUNDED 1.5 8 1 0.5 PREDICTED' % (atm1, res1, atm2, res2))
        count += 1
    if count > (seq_len * factor):
        break


# get lowest constraints and use them repulsive
"""
old_constraints.sort(key=lambda x: float(x[-1]), reverse=False)
count = 0
for constr in old_constraints:
    res1 = int(constr[0])
    res2 = int(constr[1])
    if abs(res2 - res1) >= min_dist: # and abs(res2 - res1) < rep_len * 1.5:
        atm1 = 'CB'
        if seq[res1 - 1] == 'G':
            atm1 = 'CA'
        atm2 = 'CB'
        if seq[res2 - 1] == 'G':
            atm2 = 'CA'
        score = float(constr[-1])
        score = 1 / (20 * score) # + (1/factor)
        rosetta_lines.append('AtomPair %s %d %s %d FADE -10 19 10 %.2f' % (atm1, res1, atm2, res2, round(score, 2)))
        count += 1
    if count > (seq_len * factor):
        break
"""


# write rosetta readable constraint file
#outfile_name = '.'.join(infile_name.split('.')[0:-1]) + "-" + str(factor) + '.constraints'
outfile_name = infile_name + "-" + str(factor) + '.constraints'
#outfile_name = '.'.join(infile_name.split('.')[0:-1]) + "-" + str(factor) + '.neg.constraints'
#outfile_name = '.'.join(infile_name.split('.')[0:-1]) + "-" + str(factor) + '.CA.constraints'
outfile = open(outfile_name, 'w')
for line in rosetta_lines:
    outfile.write('%s\n' % line)

#!/usr/bin/env python

import joblib, sys, os, random

infile = sys.argv[1]
name = 'pconsc26'
strict = False

if len(sys.argv) > 2:
    fraction = float(sys.argv[2])
else:
    fraction = 1

if fraction < 0:
    sys.stderr.write('Cannot use less than 0% of trees!\n')
    sys.exit()
elif fraction > 1:
    sys.stderr.write('I appreciate your enthusiasm, but I can give you only 100% of my capabilities. So I will.\n')

def predict(X, forest):
        probability = []
        for t in range(len(forest)):
            if random.random() > fraction:
                continue
        tree = forest[t]
        while len(tree) > 2:
            if X[tree[0][0]] <= tree[0][1]:
                tree = tree[1]
            else:
                tree = tree[2]
        if strict:
            if tree[1] < tree[0]:
                probability.append(0.)
            else:
                probability.append(1.)
        else:
            probability.append(tree[1]/float(tree[0] + tree[1]))
        return sum(probability)/len(probability)


if os.path.exists(infile[:infile.find('training')] +name+ '.l0'):
    sys.stderr.write('already done\n')
    sys.stderr.flush()
#   sys.exit(0)

sys.stderr.write('Forest...')
sys.stderr.flush()
forest = joblib.load('forests/layer0.dat')
X = []
Y = []
o = []
if infile.find('training') < 0:
    sys.exit(0)

of = open(infile[:infile.find('training')] + name+ '.l0', 'w')
sys.stderr.write('reading...')
sys.stderr.flush()
for l1 in open(infile):
    pair = l1[:l1.find(')')+1]
    data = l1[l1.find(')')+1:]
    data = data.split()
    addme = []
    for i in range(len(data)):
        try:
            addme.append(float(data[i]))
        except:
            pass
    pair = eval(pair)
    of.write('{:d} {:d} {:6.4f}\n'.format(pair[0], pair[1], predict(addme, forest) ) )
of.close()  


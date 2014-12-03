#!/usr/bin/env python

import joblib, sys, os, random
import numpy as np

count = 0

layer = int(sys.argv[1])
infile = sys.argv[2]
name = 'pconsc26'

strict = False


if len(sys.argv) > 3:
        fraction = float(sys.argv[3])
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

if os.path.exists(infile[:infile.rfind('training')] + name + '.l{:d}'.format(layer)):
	sys.stderr.write('already done\n')
	sys.stderr.flush()
#	sys.exit(0)

try:
	forest = joblib.load('forests/layer{:d}.dat'.format(layer))
except:
	sys.stderr.write('missing trained forest!\n')
	sys.stderr.flush()
	sys.exit(0)
X = []
Y = []
if not os.path.exists(infile[:infile.rfind('training')] + name + '.l{:d}'.format(layer-1)):
	sys.stderr.write('previous layer missing!\n')
	sys.stderr.flush()
	sys.exit(0)
previouslayer = {}
sys.stderr.write('previous layer reading...')
sys.stderr.flush()
for l1 in open(infile[:infile.rfind('training')] + name + '.l{:d}'.format(layer-1)):
	x = l1.split()
	previouslayer[(int(x[0]), int(x[1]))] = float(x[-1])
sys.stderr.write('reading initial data and predicting...')
sys.stderr.flush()
f = open(infile[:infile.rfind('training')] + name+ '.l{:d}'.format(layer), 'w')
for l1 in open(infile):
	pair = l1[:l1.find(')')+1]
	data = l1[l1.find(')') + 1:].split()

	q = []
	for qqq in data:
		q.append(eval(qqq))

	p = eval(pair)
	for i in range(-5,6):
		for j in range(-5, 6):
			try:
				q.append(previouslayer[(p[0] + i, p[1] + j)])
			except Exception as e:
				q.append(-3)
	f.write('{:d} {:d} {:6.4f}\n'.format(p[0], p[1], predict(q, forest) ) )
f.close()       
sys.stderr.write('\n')

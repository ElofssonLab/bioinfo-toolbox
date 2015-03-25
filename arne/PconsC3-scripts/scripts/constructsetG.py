#!/usr/bin/env python

# from localconfig import *
import sys
import os

if len(sys.argv) != 38:
	print 'Usage: ' + sys.argv[0] + ' <output files>'
	print 'Output files need to come in *order*!'
	print 'That is:'
	print '8x '
	print '    gplmDCA'
	print '    GaussDCA'
	print 'ML contact predictions'
	print '8x '
	print '	   Alignment stats'
	print '    Alignment'
	print 'NetSurf RSA', len(sys.argv)
	print 'SS file'
	print 'Contact file'
	print 'Outfile'
	sys.exit(1)

def parsePSIPRED(f):
	SSdict = {}
	try:
		x = open(f).read().split('\n')
	except:
		return SSdict
	for l in x:
		y = l.split()
		if len(y) != 6:
			continue
		i = int(y[0])
		SSdict[i] = [float(y[3]), float(y[4]), float(y[5])]
	return SSdict
		
def parseNetSurfP(f):
	netSurfdict = {}
	for l in open(f).readlines():
		al = []
		x = l.split()
		if l.find('#') == 0:
			continue
		if l[0] not in ['B', 'E']:
			y = ['E']
			y.extend(x)
			x = y
		for y in [4,6,7, 8, 9]:
			al.append(float(x[y]) )
		netSurfdict[ int(x[3] )] = al
	return netSurfdict

def predict(X, forest):
	probability = []
	for t in range(len(forest)):
		tree = forest[t]
		while len(tree) > 2:
			if X[tree[0][0]] <= tree[0][1]:
				tree = tree[1]
			else:
				tree = tree[2]
		probability.append(tree[1]/float(tree[0] + tree[1]))
	return sum(probability)/len(probability)



def parsePSSM(alignment):
	pssm = {}
	one2number = 'ARNDCEQGHILKMFPSTWYV-'
	bi = [ 0.0825, 0.0553, 0.0406, 0.0545, 0.0137, 0.0393, 0.0675, 0.0707, 0.0227, 0.0595, 0.0966, 0.0584, 0.0242, 0.0386, 0.0470, 0.0657, 0.0534, 0.0108, 0.0292, 0.0687 ]
	b = {}
	for i in one2number[:-1]:
		b[i] = bi[one2number.find(i)] 

	freqs = {}
	seqcount = 1.
	gapcount = 0
	coverage = []
	for l in open(alignment):
		if l.find('>') > -1:
			continue
		x = l.strip()
		if len(x) < 3:
			continue
		seqcount += 1
		coverage.append( (len(x) - x.count('-'))/float(len(x)))
		for i in range(len(x)):
			try:
				freqs[i][x[i]] += 1 
			except:
				try:
					freqs[i][x[i]] = 1 
				except:
					freqs[i] = {} 
					freqs[i][x[i]] = 1 
			if x[i] == '-':
				gapcount += 1

	b['-'] = gapcount/(seqcount * len(freqs.keys()))

	entropy = []
	for i in sorted(freqs.keys()):
		q = []
		for l in one2number:
			try:
				q.append(np.log( freqs[i][l]/(b[l] * seqcount)))
				q.append(freqs[i][l]/(b[l] * seqcount) * np.log( freqs[i][l]/(b[l] * seqcount)))
				entropy.append(freqs[i][l]/(b[l] * seqcount) * np.log( freqs[i][l]/(b[l] * seqcount)))
			except Exception as e:
				q.append(np.log( 0.1/(b[l] * seqcount)))
				q.append(0)
				entropy.append(0)
		pssm[i+1] = q
	return (pssm, np.mean(entropy), [np.min(coverage), np.max(coverage), np.mean(coverage), np.median(coverage)])

def parseStats(f):
	stats = []
	ff = open(f).readlines()
	if len(ff) != 6:
		sys.stderr.write(f + ' has incorrect format!\n')
		return [-1, -1, -1, -1, -1, -1]
	for l in ff:
		stats.append(float(l.split()[-1]))
	return stats

def parseContacts(f):
	contacts = set()
	for l in open(f):
		x = l.split()
		if len(x) != 3:
			sys.stderr.write('Incorrect format for ' + f)
			sys.exit(1)
		if float(x[-1]) < 8:
			contacts.add( (int(x[0]), int(x[1])) )
	return contacts
files = sys.argv[1:]

selected = set()
contacts = {}
X = []
Y = []
maxres = -1
acceptable = []
outfile = files[36]
if os.path.exists(outfile):
#	pass
	sys.exit(0)
import numpy as np

sys.stderr.write('Doing ' + outfile + '\n')
accessibility = parseNetSurfP(files[33])
SSdict = parsePSIPRED(files[34])
stats = []
pssm = []
entropy = []
coverage = []

for i in range(8):
	stats2  = parseStats(files[17 +i*2])
	pssm2 = parsePSSM(files[18 + i*2])
	stats.append(stats2)
	pssm.append(pssm2[0])
	entropy.append(pssm2[1])
	coverage.append(pssm2[2])

observedcontacts = parseContacts(files[35])

for index in range(17):
	contacts[index] = {}
	d = files[index]
	r = []
	if not os.path.exists(d):
		sys.stderr.write(d + ' does not exist!\n')
		sys.exit(0)
	infile = open(d).readlines()
	c = -1
	for m in infile:
		if d.find('gdca') > -1:
			x = m.split()
			divisor = 1
		elif d.find('plmdca') > -1:
			x = m.split(',')
			divisor = 1
			if len(x) != 3:
				print d + ' has wrong format!'
				sys.exit(1)
		else:
			x = m.split()
			divisor = 1
			if len(x) < 3 or x[2] != '0' or x[3] != '8':
				continue
			c = -1
		if len(x) < 3:
			continue
		aa1 = int(x[0])
		aa2 = int(x[1])
		if aa1 > maxres:
			maxres = aa1
		if aa2 > maxres:
			maxres = aa2	
#		if abs(aa1 - aa2) < 5:
#			continue
		if x[c].find('nan') > -1:
			score = -3
		else:
			score = float(x[c])
		contacts[index][(aa1, aa2)] = score/divisor
		selected.add( (aa1, aa2) )

selected2 = set()
clist = []
for c in contacts[0].keys():
	q = [ c ]
	for i in contacts.keys():
		try:
			q.append( contacts[i][c] )
		except:
			q.append( -3 )
	clist.append(q)


for i in contacts.keys():
	clist.sort(key = lambda x: -x[i+1])
	counter = -1
	c = 0
	while counter < maxres:
		j = clist[c]
		selected2.add(j[0])
		c+=1
		if abs(j[0][0] - j[0][1]) > 4:
			counter += 1

#clist.sort(key = lambda x: -1 * sum(x[1:]))
#for j in clist[:int(len(clist)*0.9)]:

# print len(selected), len(clist), float(len(selected))/len(clist)

f = open(outfile.replace('training', 'selected'), 'w')
f.write(str(selected2))
f.close()

f = open(outfile, 'w')
maxscores = []
meantop = []
stdtop = []
for index in range(17):
        maxscores.append(max(contacts[index].values()))
        q = []
        for s in list(selected2):
                try:
                        q.append(contacts[index][s])
                except:
                        pass
        meantop.append(np.mean(q))
        stdtop.append(np.std(q))

for s in sorted(list(selected)):
        q = []
        q.append(abs(s[0]-s[1]))


	for i in range(8):
		q.extend(stats[i])
		q.append(entropy[i])
		q.extend(coverage[i])
		q.extend(pssm[i][s[0]])
		q.extend(pssm[i][s[1]])

        for index in range(17):
		q.append(maxscores[index])

	for i in range(-5, 6):
                for j in range(-5, 6):
                        for index in range(17):
                                try:
                                        q.append(contacts[index][(s[0]+i, s[1]+j)])
                                        q.append((contacts[index][(s[0]+i, s[1]+j)] - meantop[index])/stdtop[index])
                                except:
                                        q.append(0)
                                        q.append(0)

        for i in range(-4, 5):
                try:
                        q.extend(SSdict[s[0]] )
                except:
                        q.extend([0,0,0])

        for i in range(-4, 5):
                try:
                        q.extend(SSdict[s[1]] )
                except:
                        q.extend([0,0,0])
	
	for i in range(-4, 5):
		try:
			q.extend(accessibility[s[0]+i] )
		except:
			q.extend([0,0,0,0,0])
#			q.extend([-1])
	
	for i in range(-4, 5):
		try:
			q.extend(accessibility[s[1]+i] )
		except:
			q.extend([0,0,0,0,0])
#			q.extend([-1])


	if s in observedcontacts:
		r = 1
	else:
		r = -1
	if True:
		X.append(q)
		Y.append(r)
#		sys.stderr.write(str(len(q)) + '\n')
		f.write(str(s))
		for qq in q:
			f.write(' ' + str(qq))
		f.write('\n')

f.close()

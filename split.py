import sys

fname = sys.argv[1]
pattern = sys.argv[2]
outprefix = sys.argv[3]

n = 0

with open(fname) as f:
    outf = open(outprefix + str(n),'w')
    for l in f:
        if pattern in l:
            n += 1
            outf.close()
            outf = open(outprefix + str(n),'w')
        outf.write(l)
    outf.close()


import sys

if len(sys.argv) < 4:
    print ("Correct usage: script.py structure1.pdb structure2.pdb output_structure.pdb")
    sys.exit()

def rewrite(infile, chain, outfile, new):
    prev_orig = ''
    with open(infile,'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                if prev_orig != '' and line[22:27] != prev_orig: new+=1
                renumb = ' '*(4-len(str(new)))+str(new)+' '
                outline = line[:21]+chain+renumb+line[27:].rstrip()
                outfile.write(outline+'\n')

            prev_orig = line[22:27]
        outfile.write('TER\n')

    return new

outfile = open(sys.argv[3],'w')

new = 1
new = rewrite(sys.argv[1], 'A', outfile, new)
new = rewrite(sys.argv[2], 'B', outfile, new+1)

outfile.write('END')
outfile.close() 
        

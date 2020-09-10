#!/usr/bin/env python3
import sys

if len(sys.argv) < 4:
    print ("Correct usage: script.py structure1.pdb structure2.pdb output_structure.pdb")
    sys.exit()


    
def rewrite(infile, chain, outfile, new):
    prev_orig = ''
    three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M'}
    
    with open(infile,'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                #print (line[22:27],":",prev_orig,":",new)
                if prev_orig != '' and line[22:27] != prev_orig:
                    new+=1
                    #prev_orig = line[22:27]
                renumb = ' '*(4-len(str(new)))+str(new)+' '
                outline = line[:21]+chain+renumb+line[27:].rstrip()
                outfile.write(outline+'\n')

            #print (line)
            prev_orig = line[22:27]
            #print ("org:",prev_orig,":",line[22:27])
        outfile.write('TER\n')

    return new

outfile = open(sys.argv[3],'w')

new = 1
new = rewrite(sys.argv[1], 'A', outfile, new)
new = rewrite(sys.argv[2], 'B', outfile, new+1)

outfile.write('END')
outfile.close() 
        

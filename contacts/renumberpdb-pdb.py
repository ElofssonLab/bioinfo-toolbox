import sys
import argparse
from math import *

# on UPPMAX only
#sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

import Bio.PDB
from Bio import pairwise2

import numpy as np

# from os.path import expanduser
# home = expanduser("~")
# sys.path.append(home + '/git/bioinfo-toolbox/')


from parsing import parse_fasta
from parsing import parse_pdb


def realign(native_filename, pdb_filename,chain,nativechains, outfilename=''):
    maxscore=0
    ### get sequence
    #seq = list(parse_fasta.read_fasta(open(fasta_filename, 'r')).values())[0][0]
    atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
    pdbfile=open(pdb_filename, 'r')
    for c in nativechains:
        seq=parse_pdb.get_atom_seq(open(native_filename, 'r'), c)
        score=pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1,score_only=True)
        if (score> maxscore):
            score=maxscore
            NC=c

    seq=parse_pdb.get_atom_seq(open(native_filename, 'r'), c)
    ref_len = len(seq)
    align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
    atom_seq_ali = align[-1][0]
    seq_ali = align[-1][1]
    #print (atom_seq_ali,seq_ali)

    res_i=-9999
    resno={}
    i=0
    atompos=0
    seqpos=0
    maxlen=len(atom_seq_ali)
    for i in range(0,maxlen):
        if atom_seq_ali[i] == "-":
            seqpos+=1
        elif seq_ali[i] =="-":
            atompos+=1
            resno[atompos]=-9999
        else:
            atompos+=1
            seqpos+=1
            resno[atompos]=seqpos
    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)
    i=0
    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_pdb.parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain  and chain != '*':
            continue
        if atm_record['res_no'] != res_i:
            i+=1
            res_i = atm_record['res_no']
        atm_record['res_no']=resno[i]
        #print (atm_record)
        if resno[i]>0:
            parse_pdb.write_pdb_atm_record(atm_record)
        #res_dict[res_i].append(np.array(atm))
        

    #pdbfile.close()

    return(NC)
  
    



if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('native')#, required=True)
    p.add_argument('pdb')
    p.add_argument('-o', '--outfile', default='')
    #p.add_argument('-c', '--chain',default='*')
    #p.add_argument('-n', '--nativechain',default='*')
    args = vars(p.parse_args(sys.argv[1:]))


    
    #if len(open(args['pdb']).readline().split(' ')) != 3:

    #print (args['pdb'],args['native'])
    modelchains=parse_pdb.get_all_chains(args['pdb'])
    nativechains=parse_pdb.get_all_chains(args['native'])

    #print ("Chains",modelchains,nativechains)
    for chain in modelchains:
        #print ("TEST",chain)
        NC=realign(args['native'],  args['pdb'],chain,nativechains,outfilename=args['outfile'])
        nativechains.remove(NC)

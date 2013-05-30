import sys
import os
import parse_fasta
import shutil
import errno
import stat

pdb_path = '/bubo/home/h9/mircomic/glob/databases/PDB/current_release/pdb_seqres.txt'
num_rep = -1

pdb_file = open(pdb_path, 'r')
pdb_seq_dict = parse_fasta.read_fasta_pdb(pdb_file)
pdb_file.close()

flist_string = ''
query_file = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

for line in query_file:
    line_arr = line.split('\t')
    
    #acc = line_arr[0].split('_')[0]
    acc = line_arr[0]
    print acc

    if acc[0] == '#':
        continue
    
    chain = line_arr[0].split('_')[1]

    #pos = (int(line_arr[1].split('-')[0]), int(line_arr[1].split('-')[1]))
    pos = eval(line_arr[1])
    if num_rep != -1 and len(pos) < num_rep:
        continue

    name = line_arr[2]

    seq = pdb_seq_dict[acc][0]

    if num_rep != -1:
        seq_start = pos[0][0]
        seq_end = pos[num_rep - 1][1]
        query_seq = seq[seq_start-1:seq_end]
    else:
        seq_start = 1
        seq_end = len(seq)
        query_seq = seq

    fasta_string = '>%s/%s/%d-%d\n%s\n' % (name, acc, seq_start, seq_end, query_seq)
    #print(fasta_string)

    outfile.write(fasta_string)


query_file.close()
outfile.close()

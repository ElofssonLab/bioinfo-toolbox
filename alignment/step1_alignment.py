import sys
import os
import parse_fasta

pdb_path = '/bubo/home/h9/mircomic/glob/databases/PDB/current_release/pdb_seqres.txt'
hhblits_path = './loop_hhblits3/'
jhmmer_path = './loop_jackhmmer/'

num_rep = 2

pdb_file = open(pdb_path, 'r')
pdb_seq_dict = parse_fasta.read_fasta_pdb(pdb_file)
pdb_file.close()

flist_string = ''
query_file = open(sys.argv[1], 'r')

for line in query_file:
    line_arr = line.split('\t')
    
    #acc = line_arr[0].split('_')[0]
    acc = line_arr[0]
    chain = line_arr[0].split('_')[1]
    #pos = (int(line_arr[1].split('-')[0]), int(line_arr[1].split('-')[1]))
    pos = eval(line_arr[1])
    name = line_arr[2]

    seq = pdb_seq_dict[acc][0]
    seq_start = pos[0][0]
    seq_end = pos[num_rep - 1][1]
    query_seq = seq[seq_start-1:seq_end]

    flist_string += '%s_r%d.fa\n' % (acc, num_rep)
    fasta_string = '>%s/%d-%d\n%s\n' % (acc, seq_start, seq_end, query_seq)
    #print(fasta_string)

    input_hhblits_file = open('%s/input/%s_r%d.fa' % (hhblits_path, acc, num_rep), 'w')
    input_hhblits_file.write(fasta_string)
    input_hhblits_file.close()

    input_jhmmer_file = open('%s/input/%s_r%d.fa' % (jhmmer_path, acc, num_rep), 'w')
    input_jhmmer_file.write(fasta_string)
    input_jhmmer_file.close()

print flist_string
flist_hhblits_file = open('%sflist.txt' % hhblits_path, 'w')
flist_hhblits_file.write(flist_string)
flist_hhblits_file.close()

#flist_jhmmer_file = open('%sflist.txt' % jhmmer_path, 'w')
#flist_jhmmer_file.write(flist_string)
#flist_jhmmer_file.close()

os.chdir(hhblits_path)
hhblits_cmd = 'python wrap_run_hhblits_and_ECs.py flist.txt'
os.system(hhblits_cmd)

#os.chdir('..')
#os.chdir(jhmmer_path)
#jhmmer_cmd = 'python wrap_run_jackhammer_and_ECs.py flist.txt'
#os.system(jhmmer_cmd)
os.chdir('..')

import sys
import os
import parse_fasta

num_rep = 2
evals = [3, 5, 10, 15, 20, 30, 40]
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
    query_seq = seq[pos[0][0]-1:pos[1][1]]

    flist_string += '%s_r2.fa\n' % acc
    fasta_string = '>%s/%d-%d\n%s\n' % (acc, pos[0][0], pos[1][1], query_seq)
    #print(fasta_string)

    input_hhblits_file = open('%s/input/%s_r2.fa' % (hhblits_path, acc), 'w')
    input_hhblits_file.write(fasta_string)
    input_hhblits_file.close()

    input_jhmmer_file = open('%s/input/%s_r2.fa' % (jhmmer_path, acc), 'w')
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

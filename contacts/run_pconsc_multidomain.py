import sys
import os
import parse_fasta
import shutil
import errno
import stat

pdb_path = '/bubo/home/h9/mircomic/glob/databases/PDB/current_release/pdb_seqres.txt'
main_path = '/bubo/home/h9/mircomic/glob/project_multidomain/pconsc_predictions'

pdb_file = open(pdb_path, 'r')
pdb_seq_dict = parse_fasta.read_fasta_pdb(pdb_file)
pdb_file.close()

flist_string = ''
query_file = open(sys.argv[1], 'r')

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

    seq = pdb_seq_dict[acc][0]
    seq_start = pos[0][0]
    seq_end = pos[-1][1]
    query_seq = seq[seq_start-1:seq_end]

    fasta_string = '>%s/%d-%d\n%s\n' % (acc, seq_start, seq_end, query_seq)
    #print(fasta_string)

    try:
        os.mkdir('%s/%s' % (main_path, acc))
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir('%s/%s' % (main_path, acc)):
            pass
        else: raise

    os.chdir('%s/%s' % (main_path, acc))

    ### run pconsc on given sequence (submit job)
    """
    seq_file = open('%s.fa' % acc, 'w')
    seq_file.write(fasta_string)
    seq_file.close()
    shutil.copyfile('../run_pconsc.sh', './run_pconsc.sh')
    st = os.stat('./run_pconsc.sh')
    os.chmod('./run_pconsc.sh', st.st_mode | stat.S_IEXEC)
    os.system('sbatch ./run_pconsc.sh %s.fa' % acc)
    """

    ### run psipred on given sequence (submit job)
    """
    os.system('sbatch /bubo/home/h9/mircomic/glob/psipred/run_psipred.sh %s.fa &' % acc)
    """

    ### remove fasta alignment files
    """
    filelist = [ f for f in os.listdir(".") if f.endswith(".fas") ]
    for f in filelist:
        os.remove(f)
        print 'removing %s' % f
    """


    ### get reference PDB structures for contact map plotting
    """
    print acc.split('_')[0]
    print acc.split('_')[0][1:3]

    if not os.path.exists('pdb%s.ent.gz' % acc.split('_')[0]) and not os.path.exists('pdb%s.ent' % acc.split('_')[0]):
        os.system('wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz' % (acc.split('_')[0][1:3], acc.split('_')[0]))
    if not os.path.exists('pdb%s.ent' % acc.split('_')[0]):
        os.system('gunzip pdb%s.ent.gz' % acc.split('_')[0])
    """

    ### plot contact map of pconsc output
    """
    try:
        shutil.move('pconsc.slurm.out', '%s.out' % acc)
    except IOError as exc:
        if os.path.exists('%s.out' % acc):
            pass
        else: raise
    os.system('python ../../plot_contact_map.py %s.out %s.horiz pdb%s.ent 1' % (acc, acc, acc.split('_')[0]))
    shutil.move('%s.out_ContactMap.pdf' % acc, '../')
    #shutil.move('%s.out_PPVs.pdf' % acc, '../')
    """




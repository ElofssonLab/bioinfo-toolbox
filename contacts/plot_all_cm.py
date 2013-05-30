import os

for root, dirs, filenames in os.walk('pconsc_predictions'):
    for f in filenames:
        #if f.endswith('plmdca') or f.endswith('psicov'):
        if f.endswith('contacts') and len(f.split('.')) == 2:
            acc_chain = f.split('.')[0]
            acc = acc_chain.split('_')[0]
            cm_file = os.path.join(root, f)
            fasta_file = os.path.join(root, '%s.fa' % acc_chain)
            pdb_file = os.path.join(root, 'pdb%s.ent' % acc)
            print acc_chain
            print cm_file
            os.system('python /bubo/home/h9/mircomic/glob/project_multidomain/plot_contact_map.py %s %s %s 1000' % (cm_file, fasta_file, pdb_file))

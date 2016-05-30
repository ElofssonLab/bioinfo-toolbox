import sys
from collections import defaultdict


# Parser for CATH list files (clf)
#
# CATH IDs: 
#  - superfamily = X.X.X.X
#  - domain = '4-letter-PDB' + 'chain' + 'domain', eg: 1xxxA01


# returns a dictionary {superfamily_id : domain_id_list}
def get_sfam_dom_dict(clf):

    cath_dict = defaultdict(list)

    for line in clf:
        if line.startswith('#'):
            continue
        line_arr = line.split()
        
        cath_dom = line_arr[0]
        cath_sfam = '.'.join([line_arr[1], line_arr[2], line_arr[3], line_arr[4]])

        cath_dict[cath_sfam].append(cath_dom)

    return cath_dict


# returns a dictionary {domain_id : superfamily_id_list}
def get_dom_sfam_dict(clf):

    cath_dict = defaultdict(list)

    for line in clf:
        if line.startswith('#'):
            continue
        line_arr = line.split()
        
        cath_dom = line_arr[0]
        cath_sfam = '.'.join([line_arr[1], line_arr[2], line_arr[3], line_arr[4]])

        cath_dict[cath_dom].append(cath_sfam)

    return cath_dict


# returns a dictionary {pdb_id : superfamily_id_list}
def get_pdb_sfam_dict(clf):

    cath_dict = defaultdict(list)

    for line in clf:
        if line.startswith('#'):
            continue
        line_arr = line.split()
        
        pdb_chain_id = line_arr[0][:4]
        cath_sfam = '.'.join([line_arr[1], line_arr[2], line_arr[3], line_arr[4]])

        cath_dict[pdb_chain_id].append(cath_sfam)

    return cath_dict



if __name__ == '__main__':

    test_fam = get_sfam_dom_dict(open(sys.argv[1], 'r'))
    test_dom = get_dom_sfam_dict(open(sys.argv[1], 'r'))
    test_pdb = get_pdb_sfam_dict(open(sys.argv[1], 'r'))
    
    for acc, sfam in test_pdb.iteritems():
        print acc, sfam
    #print len(test_fam)
    #print len(test_dom)
    #print len(test_pdb)

    #print test_fam['1.10.8.10']
    #print test_dom['3fe3A03']
    #print test_pdb['3fe3']

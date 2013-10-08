import sys

import parse_cath_clf

if __name__ == '__main__':

    cath_dict = parse_cath_clf.get_pdb_sfam_dict(open('/bubo/home/h9/mircomic/glob/databases/CATH/CathDomainList.S100.v3.5.0', 'r'))
    
    for id in open('IDs_all.txt', 'r'):
        id = id.strip()
        cath_id = cath_dict[id.lower()]
        if(len(cath_id) > 0):
            print id + ':A ' + ' '.join(cath_id[0].split('.'))
        else:
            print id + ':A '

        #print id
        #print cath_id
        #print id + ':A\t' + '\t'.join(cath_id.split('.'))

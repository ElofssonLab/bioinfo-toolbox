#!/bin/python
import sys
sys.path.append('/home/mircomic/toolbox')
from parsing import parse_vall
from parsing import parse_cath_clf
from time import clock


def filter(cath_dom_list, frag_dict):

    count = 0

    for dom_id in cath_dom_list:
        pdb_id = dom_id[:4]
        if pdb_id in frag_dict:
            count += 1
            del frag_dict[pdb_id]

    return frag_dict, count



def main(vallfile_name, clf_name, pdb_id):

    pdb_sfam_dict = parse_cath_clf.get_pdb_sfam_dict(open(clf_name, 'r'))
    sfam_dom_dict = parse_cath_clf.get_sfam_dom_dict(open(clf_name, 'r'))

    sfam_list = pdb_sfam_dict[pdb_id]
    dom_list = []

    for sfam in sfam_list:
        dom_list += sfam_dom_dict[sfam]

    print dom_list
    
    print 'Parsing started...'
    t0 = clock()
    vallfile_list = parse_vall.read(open(vallfile_name, 'r'))
    t1 = clock()
    print 'Parsing ended in %ds.\n' % (t1 - t0)

    frag_dict = vallfile_list[2]
    
    print 'Running filter...'
    t0 = clock()
    new_frag_dict, count = filter(dom_list, frag_dict)
    t1 = clock()
    print 'Filtered in %ds.\n' % (t1 - t0)
    vallfile_list[2] = new_frag_dict

    return vallfile_list, count


if __name__ == '__main__':

    new_vallfile_list, count = main(sys.argv[1], sys.argv[2], sys.argv[3])
    print 'Removed %d homologous proteins from vall file.' % count

    print 'Writing started...'
    t0 = clock()
    parse_vall.write(new_vallfile_list, open('%s' % sys.argv[1], 'w'))
    t1 = clock()
    print 'Writing ended in %ds.\n' % (t1 - t0)


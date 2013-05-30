import sys
import parse_fragments
import parse_cath_clf


def filter(cath_dom_list, fragfile_list, frag_len):

    new_fragfile_list = []
    counts = []

    for (pos_header, frag_dict) in fragfile_list:
        old_num_frag = int(pos_header.split()[-1])
        new_num_frag = int(pos_header.split()[-1])
        for dom_id in cath_dom_list:
            pdb_id = dom_id[:4]
            if frag_dict.has_key(pdb_id):
                new_num_frag -= len(frag_dict[pdb_id]) / frag_len
                del frag_dict[pdb_id]
        counts.append(new_num_frag)

        if new_num_frag < 100 and new_num_frag >= 10:
            new_num_frag_str = ' %d' % new_num_frag
            new_header = pos_header.replace(str(old_num_frag),new_num_frag_str)
        if new_num_frag < 10:
            new_num_frag_str = '  %d' % new_num_frag
            new_header = pos_header.replace(str(old_num_frag), new_num_frag_str)
        if new_num_frag <= 0:
            sys.stderr.write('All fragments for position %s have been filtered.')
        else:
            new_header = pos_header.replace(str(old_num_frag), str(new_num_frag))
        new_fragfile_list.append((new_header, frag_dict))

    return new_fragfile_list, counts



def main(fragfile_name, clf_name, pdb_id, chain):

    pdb_sfam_dict = parse_cath_clf.get_pdb_sfam_dict(open(clf_name, 'r'))
    sfam_dom_dict = parse_cath_clf.get_sfam_dom_dict(open(clf_name, 'r'))

    pdb_chain_id = pdb_id + chain
    #sfam_list = pdb_sfam_dict[pdb_chain_id]
    sfam_list = pdb_sfam_dict[pdb_id]
    dom_list = []

    #print sfam_list

    for sfam in sfam_list:
        dom_list += sfam_dom_dict[sfam]

    #print dom_list

    fragfile_list = parse_fragments.read(open(fragfile_name, 'r'))
    frag_len = parse_fragments.get_frag_len(open(fragfile_name, 'r'))
    
    new_fragfile_list, counts = filter(dom_list, fragfile_list, frag_len)

    return new_fragfile_list, counts


if __name__ == '__main__':

    new_fragfile_list, counts = main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    
    #print len(counts)
    #print counts
    print 'Minimum number of fragments at any position: %d' % min(counts)
    
    frag_len = parse_fragments.get_frag_len(open(sys.argv[1], 'r'))
    parse_fragments.write(new_fragfile_list, open('%s' % sys.argv[1], 'w'), frag_len)


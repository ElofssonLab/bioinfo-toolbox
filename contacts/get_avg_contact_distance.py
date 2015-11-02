import sys
import math

import parse_contacts
import parse_fasta


def get_min_dist(contact, contact_lst):

    min_dist = float('inf')

    x0 = contact[1]
    y0 = contact[2]

    for curr_contact in contact_lst:
        x1 = curr_contact[1]
        y1 = curr_contact[2]
        if x0 == x1 and y0 == y1:
            continue
        dist = math.sqrt(math.pow(x1-x0, 2) + math.pow(y1-y0, 2))
        if dist < min_dist:
            min_dist = dist
            min_x = x1
            min_y = y1

    #print '(%d,%d) to (%d,%d): %d' % (x0, y0, min_x, min_y, min_dist)
    
    return min_dist


def main(seqfile, contactfile, sep=' '):
    
    seq = parse_fasta.read_fasta(seqfile).values()[0][0]
    seq_len = len(seq)

    all_contacts = parse_contacts.parse(contactfile, sep)
    top_contacts = []
    for i, contact in enumerate(all_contacts):
        if abs(contact[1] - contact[2]) < 10:
            continue
        top_contacts.append(contact)
        if i >= seq_len:
            break
    min_dist_lst = []
    for contact in top_contacts:
        min_dist = get_min_dist(contact, top_contacts)
        min_dist_lst.append(min_dist)
    #return sum(min_dist_lst) / float(len(min_dist_lst))
    return min_dist_lst



if __name__ == '__main__':
    
    idfile = open(sys.argv[1], 'r')
    dist_dict = {}
    for line in idfile:
        id = line.strip().split()
        if len(id) == 1:
            id = id[0]
            seqfile = open('data/%s:A/sequence.fa' % id, 'r')
            contactfile = open('data/%s:A/pconse.results' % id, 'r')
            dist_dict[id] = main(seqfile, contactfile)
        else:
            seqfile = open('pfam/data/pconsc_predictions/%s/%s.fa' % (id[0], id[1]), 'r')
            contactfile = open('pfam/data/pconsc_predictions/%s/%s.out' % (id[0], id[1]), 'r')
            dist_dict[id[0]] = main(seqfile, contactfile)
        seqfile.close()
        contactfile.close()
        #contactfile = open('data/%s:A/mar23-hh0hh4hh10hhb40jh0jh4jh10jhm40-psicovplmdca.results' % id, 'r')


    idfile.close()

    id_dist_lst = sorted(dist_dict.iteritems(), key=lambda(x):x[0])
    #print id_dist_lst
    outfile = open(sys.argv[2], 'w')
    for (id, dist_lst) in id_dist_lst:
        outfile.write('%s\t%s\n' % (id, '\t'.join(map(str, dist_lst))))
    outfile.close()




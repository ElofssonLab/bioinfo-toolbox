import sys

import parse_fasta
import parse_contacts

if __name__ == '__main__':

    idfile = open(sys.argv[1], 'r')
    score_dict = {}

    for line in idfile:
        id = line.strip()

        seqfile = open('data/%s:A/sequence.fa' % id, 'r')
        #contactfile = open('data/%s:A/mar23-hh0hh4hh10hhb40jh0jh4jh10jhm40-psicovplmdca.results' % id, 'r')
        contactfile = open('data/%s:A/pconse.results' % id, 'r')

        seq = parse_fasta.read_fasta(seqfile).values()[0][0]

        ref_len = len(seq)
        contacts = parse_contacts.parse(contactfile, ' ')
        if contacts == []:
            continue
        score_lst = []
        i = 0
        for contact in contacts:
            score_lst.append(contact[0])
            i += 1
            if i >= ref_len * 1.0:
                break
        avg_score = sum(score_lst) / float(len(score_lst))
        score_dict[id] = avg_score

    idfile.close()
    
    id_score_lst = sorted(score_dict.iteritems(), key=lambda(x):x[0])
    outfile = open(sys.argv[2], 'w')
    for (id, score) in id_score_lst:
        outfile.write('%s\t%s\n' % (id, score))
    outfile.close()

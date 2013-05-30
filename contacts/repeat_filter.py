import sys

import parse_contacts
import parse_fasta
import parse_hmm
import dotter


def filter_contacts(contact_filename, fasta_filename, sep=','):  
    
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ss = seq
    ref_len = len(seq)

    contacts = parse_contacts.parse(open(contact_filename, 'r'), sep)
    profile = parse_hmm.read_hmm(open('%s.hmm' % '.'.join(contact_filename.split('.')[:-1])))
    dot_matrix = dotter.calc_dot_matrix_profile(seq, profile)

    new_contacts = []
    count = 0
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1

        pos_sim = dot_matrix[c_x, c_y]

        #new_score = max(score * (1 - pos_sim), 0)
        if pos_sim > 0.1:
            new_score = 0.0
        else:
            new_score = score
        new_contacts.append((new_score, c_x + 1, c_y + 1))

        count += 1
        #if count > ref_len * 1.0:
        #    break

    new_contacts.sort(key=lambda x: x[0], reverse=True)

    outfilename = '%snew' % contact_filename
    parse_contacts.write(new_contacts, open(outfilename, 'w'), sep)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        sys.stderr.write('Usage: python repeat_filter.py <contact_filename> <fasta_filename>\n')

    contact_filename = sys.argv[1]
    fasta_filename = sys.argv[2]

    # guessing separator of constraint file
    line = open(contact_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'

    filter_contacts(contact_filename, fasta_filename, sep)


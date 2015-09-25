#!/usr/bin/env python

"""Contact Map Filter: refines contact maps by adjusting contact scores.

Scores are either set to zero for impossible/improbable contacts, such as:
- contacts within the same secondary structural element
- contacts that are close in sequence
- cysteine contacts to multiple residues

or set to one in case of:
- top ranked cys-cys contacts2
"""

import sys


#For some reason relative paths did not work on my mac /AE
from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox/parsing')
sys.path.append(home + '/git/bioinfo-toolbox/parsing')
sys.path.append(home + '/bioinfo-toolbox/contacts')
sys.path.append(home + '/git/bioinfo-toolbox/contacts')
 
import parse_contacts 
import parse_psipred 
import parse_fasta 
import plot_contact_map


def get_ss_pos(ss_seq):

    """Get positions of secondary structure elements in sequence
    @param  ss_seq  string of secondary structure predictions
    @return [(start_pos, end_pos)]
    """

    ss_pos_lst = []
    prev_ss = ss_seq[0]
    start = 1
    end = len(ss_seq)
    in_ss_flag = False

    for i, ss in enumerate(ss_seq):
        if ss != 'C' and ss != prev_ss: # change in ss => new element
            if in_ss_flag: # no coil between secstruct elements
                end = i
                ss_pos_lst.append((start, end))
            start = i + 1
            in_ss_flag = True
        elif in_ss_flag and ss != prev_ss: # change + already in ss => element ends
            end = i
            ss_pos_lst.append((start, end))
            in_ss_flag = False
        prev_ss = ss

    if len(ss_pos_lst) == 0:
        ss_pos_lst.append(start, end)
    return ss_pos_lst


def secstruct_filter(c_lst, ss_seq):

    """Filters all contacts within the same secondary structural element
    @param  c_lst   contact list (as given by parsing/parse_contacts.py)
    @param  ss_seq  string of secondary structure predictions
    Ensures: len(c_lst) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """
   
    c_filt = []
    ss_pos_lst = get_ss_pos(ss_seq)
    #print ss_pos_lst
    for c in c_lst:
        score = c[0]
        res1 = c[1]
        res2 = c[2]
        for ss_pos in ss_pos_lst:
            start = ss_pos[0]
            end = ss_pos[1]
            if (res1 >= start and res2 <= end) and (res1 <= end and res2 >= start):
                #print '%d - %d: %f in %d - %d' % (res1, res2, score, start, end)
                score = 0.0
                break
            if (res2 >= start and res1 <= end) and (res2 <= end and res1 >= start):
                #print '%d - %d: %f in %d - %d' % (res1, res2, score, start, end)
                score = 0.0
                break
        c_filt.append((score, res1, res2))

    return c_filt


def conservation_filter(c_lst, ali):

    """Filters contacts from highly conserved residues
    @param  c_lst       contact list (as given by parsing/parse_contacts.py)
    @param  ali     multiple sequence alignment ([seq_1, ..., seq_n])
    Ensures: len(c_lst) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """

    c_filt = []
    return c_filt


def seq_filter(c_lst, i):

    """Filters all contacts that are less than i residues apart in the sequence
    @param  c_lst   contact list (as given by parsing/parse_contacts.py)
    @param  i   integer sequence distance
    Ensures: len(c_lst) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """

    c_filt = []
    return c_filt


def cysteine_filter(c_lst, seq):

    """Filters contacts from cysteines
    @param  c_lst       contact list (as given by parsing/parse_contacts.py)
    @param  seq     string of one-letter coded amino acid sequence
    Ensures: len(c_lst) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """

    c_filt = []
    return c_filt


if __name__ == "__main__":

    c_filename = sys.argv[1]
    psipred_filename = sys.argv[2]
    cfilt_filename = sys.argv[3]
    #seq_filename = sys.argv[3]

    c_lst = parse_contacts.parse(open(c_filename, 'r'))
    ss_seq = parse_psipred.horizontal(open(psipred_filename, 'r'))
    c_filt = secstruct_filter(c_lst, ss_seq)

    cfilt_file = open(cfilt_filename, 'w')
    parse_contacts.write(c_filt, cfilt_file)
    cfilt_file.close()

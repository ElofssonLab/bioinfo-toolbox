import sys

def secstruct_filter(c, ss):

    """Filters all contacts within the same secondary structural element
    @param  c   contact list (as given by parsing/parse_contacts.py)
    @param  ss  string of secondary structure predictions
    Ensures: len(c) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """
   
    c_filt = []
    return c_filt


def conservation_filter(c, ali):

    """Filters contacts from highly conserved residues
    @param  c       contact list (as given by parsing/parse_contacts.py)
    @param  ali     multiple sequence alignment ([seq_1, ..., seq_n])
    Ensures: len(c) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """

    c_filt = []
    return c_filt


def seq_filter(c, i):

    """Filters all contacts that are less than i residues apart in the sequence
    @param  c   contact list (as given by parsing/parse_contacts.py)
    @param  i   integer sequence distance
    Ensures: len(c) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """

    c_filt = []
    return c_filt


def cysteine_filter(c, seq):

    """Filters contacts from cysteines
    @param  c       contact list (as given by parsing/parse_contacts.py)
    @param  seq     string of one-letter coded amino acid sequence
    Ensures: len(c) == len(c_filt), only contact weights are changed
    @return [(score, residue a, residue b)]
    """

    c_filt = []
    return c_filt


#!/usr/bin/env python


def parse(afile, sep=' ', min_dist=4):
    
    """Parse contact file.
    @param  afile   contact file
    @param  sep     separator of contact file (default=' ')
    Ensures: Output is sorted by confidence score.
    @return [(score, residue a, residue b)]
    """

    contacts = []
    for aline in afile:
        if aline.strip() != '':
            line_arr = aline.strip().split(sep)
            if line_arr[0].startswith('E'):
                continue
            i = int(line_arr[0])
            j = int(line_arr[1])
            score = float(line_arr[-1])
            if abs(i - j) > min_dist:
                contacts.append((score, i, j))
    afile.close()

    contacts.sort(key=lambda x: x[0], reverse=True)
    return contacts


def write(contacts, outfile, sep=' '):

    """Write contact file.
    @param  contacts    contact list
    @param  outfile     output contact file
    @param  sep     separator of contact file (default=' ')
    """

    for c in contacts:
        outfile.write('%d%s%d%s%f\n' % (c[1], sep, c[2], sep, c[0]))


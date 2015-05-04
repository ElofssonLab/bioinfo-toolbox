import sys
import os
import shutil
import errno
import stat

from parsing import parse_fasta

if __name__ == '__main__':

    infile = open(sys.argv[1], 'r')
    seq_dict = parse_fasta.read_fasta(infile)

    for header, seq_lst in seq_dict.iteritems():
        acc = header.split()[0]
        outfile = open(acc + '.fa', 'w')

        seq = seq_lst[0]
        fasta_string = '>%s\n%s\n' % (header, seq)
        #print(fasta_string)

        outfile.write(fasta_string)
        outfile.close()


    infile.close()

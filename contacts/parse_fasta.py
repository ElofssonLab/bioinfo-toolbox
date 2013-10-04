import string, copy
import sys

## parses any fasta, a2m, a3m file, sequence or alignment
## ensures: key of a given query ID only contains its ID, not the full header
## output: {header: [sequence_1, sequence_2, ...]}
def read_fasta(afile, query_id=''):

    seq_dict = {}

    header = ''
    seq = ''

    for aline in afile:
        aline = aline.strip()

        # check for header
        if aline.startswith('>'):
            if header != '' and seq != '':
                if seq_dict.has_key(header):
                    seq_dict[header].append(seq)
                else:
                    seq_dict[header] = [seq]
            seq = ''
            if aline.startswith('>%s' % query_id) and query_id !='':
                header = query_id
            else:
                header = aline[1:]

        # otherwise concatenate sequence
        else:
            #aline_seq = aline.translate(None, '.-').upper()
            seq += aline

    # add last entry
    if header != '':
        if seq_dict.has_key(header):
            seq_dict[header].append(seq)
        else:
            seq_dict[header] = [seq]
    else:
        sys.stderr.write('ERROR: file empty or wrong file format')

    return seq_dict


## parses any fasta, a2m, a3m file, sequence or alignmen
## ensures: key = PDB accession
## output: {PDB-acc: [sequence_1, sequence_2, ...]}
def read_fasta_pdb(afile, query_id=''):

    seq_dict = {}

    header = ''
    seq = ''

    for aline in afile:
        aline = aline.strip()

        # check for header
        if aline.startswith('>'):
            if header != '' and seq != '':
                if seq_dict.has_key(header):
                    seq_dict[header].append(seq)
                else:
                    seq_dict[header] = [seq]
            seq = ''
            if aline.startswith('>%s' % query_id) and query_id !='':
                header = query_id
            else:
                header = aline[1:].split()[0]

        # otherwise concatenate sequence
        else:
            #aline_seq = aline.translate(None, '.-').upper()
            seq += aline

    # add last entry
    if header != '':
        if seq_dict.has_key(header):
            seq_dict[header].append(seq)
        else:
            seq_dict[header] = [seq]
    else:
        sys.stderr.write('ERROR: file empty or wrong file format')

    return seq_dict



if __name__ == "__main__":

    afile = open(sys.argv[1], 'r')
    if len(sys.argv) == 3:
        query_id = sys.argv[2]
    else:
        query_id = ''
    seq_dict = read_fasta(afile, query_id)
    afile.close()

    print seq_dict
    #print seq_dict['>sp|P20749|BCL3_HUMAN B-cell lymphoma 3 protein OS=Homo sapiens GN=BCL3 PE=1 SV=2']
    #print len(seq_dict)

from pymol import cmd

import sys

sys.path.append('/home/mirco_local/bioinfo-toolbox')
print sys.path
#from parsing import parse_contacts
#from parsing import parse_fasta
#import parse_contacts
#import parse_fasta

def parse_contacts(afile, sep=' '):
    
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
            if abs(i - j) > 4:
                contacts.append((score, i, j))
    afile.close()

    contacts.sort(key=lambda x: x[0], reverse=True)
    return contacts


def read_fasta(afile, query_id=''):

    """Parses any fasta, a2m, a3m file, sequence or alignment file.
    @param  afile       input file
    @param  query_id    ID of query sequence (default='')
    Ensures: key of a given query ID only contains its ID, not the full header
    @return {header: [sequence_1, sequence_2, ...]} 
    """

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


prot_name = sys.argv[1]
seqfile = sys.argv[2]
pdbfile = sys.argv[3]
cfile = sys.argv[4]
factor = float(sys.argv[5])

seqlen = len(read_fasta(open(seqfile, 'r')).values()[0][0])
contacts = parse_contacts(open(cfile, 'r'))


cmd.load(pdbfile)
cmd.set('dash_gap', 0.0)
cmd.set('dash_radius', 0.1)
cmd.bg_color('white')
cmd.hide('everything')
cmd.show('cartoon')
cmd.color('gray', prot_name)
cmd.set('ray_shadows', 0)

view = (\
    -0.929645956,   -0.144559741,   -0.338900387,\
     0.319600046,   -0.774050355,   -0.546532333,\
    -0.183317721,   -0.616394758,    0.765792787,\
     0.000000000,    0.000000000, -157.942276001,\
    18.551437378,    6.361377716,   -4.695641994,\
   131.235534668,  184.649002075,  -20.000000000 )

cmd.set_view(view)

count = 0

for contact in contacts:
    if count > seqlen * factor:
        break
    atm1 = '/%s//A/%s/CA' % (prot_name, contact[1])
    atm2 = '/%s//A/%s/CA' % (prot_name, contact[2])
    cmd.select('a', atm1)
    cmd.select('b', atm2)
    d_name = 'd' + str(count)
    cmd.distance(d_name, atm1, atm2)
    cmd.hide('labels', d_name)

    col_name = 'own_color' + str(count)
    score = contact[0]
    dist = cmd.distance(d_name, atm1, atm2)
    #if score > 0.5:
    if dist < 10:
        #cmd.hide('labels', d_name)
        #cmd.set_color(col_name, [0.0 + score, 1.0 - score, 0.0])
        #cmd.color(col_name, d_name)
        cmd.color("green", d_name)
    elif dist < 12:
        cmd.color("orange", d_name)
    else:
        #cmd.set_color(col_name, [0.0, 0.1 + score,1.0 - score])
        #cmd.color(col_name, d_name)
        cmd.color("red", d_name)

    count = count +1





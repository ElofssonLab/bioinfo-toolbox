import sys, os

 
# pymol environment
moddir='/opt/pymol-svn/modules'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')

from pymol import cmd


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


def draw(prot_name, seqfile, pdbfile, cfile, factor):
    seq = read_fasta(open(seqfile, 'r')).values()[0][0]
    seqlen = len(seq)
    contacts = parse_contacts(open(cfile, 'r'))

    cmd.load(pdbfile)
    cmd.set('dash_gap', 0.0)
    cmd.set('dash_radius', 0.1)
    cmd.bg_color('white')
    cmd.hide('everything')
    cmd.show('cartoon')
    cmd.color('gray', prot_name)
    
    """
    view = (\
        -0.271442711,   -0.905138493,    0.327085078,\
         0.034359235,   -0.348747194,   -0.936563492,\
         0.961805284,   -0.242983535,    0.125761271,\
         0.000000000,    0.000000000, -128.858474731,\
        17.616867065,   -0.161788940,   -4.633638382,\
       103.167495728,  154.549499512,  -20.000000000 )

    #PconsFold paper:
        -0.087676540,   -0.441159278,   -0.893107057,\
         0.482956320,    0.765313327,   -0.425448418,\
         0.871214509,   -0.468636513,    0.145965248,\
        -0.000171857,    0.001094781, -158.355941772,\
        23.238109589,   -7.513275623,    6.870838165,\
       131.660964966,  185.186004639,  -20.000000000 )
    cmd.set_view(view)
    """

    count = 0
    maxdist = 20.0

    for contact in contacts:
        if count > seqlen * factor:
            break
        if seq[contact[1]-1] == 'G':
            atm1a = '/%s//A/%s/CA' % (prot_name, contact[1])
            atm1b = atm1a
        else:
            atm1a = '/%s//A/%s/CA' % (prot_name, contact[1])
            atm1b = '/%s//A/%s/CB' % (prot_name, contact[1])
        if seq[contact[2]-1] == 'G':
            atm2a = '/%s//A/%s/CA' % (prot_name, contact[2])
            atm2b = atm2a
        else:
            atm2a = '/%s//A/%s/CA' % (prot_name, contact[2])
            atm2b = '/%s//A/%s/CB' % (prot_name, contact[2])
        cmd.select('a', atm1a)
        cmd.select('b', atm2a)
        d_name = 'd' + str(count)
        dist = cmd.distance(d_name, atm1a, atm2a)
        if dist == -1:
            ## atom pair not in structure (xtal vs. seqres)
            continue
        cmd.hide('labels', d_name)

        col_name = 'own_color' + str(count)
        score = contact[0]
        dist = cmd.distance(d_name + 'b', atm1b, atm2b)
        cmd.delete(d_name + 'b')
        #if score > 0.5:
        print '%s: %f => %f' % (d_name, dist, dist/(maxdist/2))
        if dist < maxdist/2:
            #cmd.hide('labels', d_name)
            #cmd.set_color(col_name, [0.0 + score, 1.0 - score, 0.0])
            cmd.set_color(col_name, [0.0 + (dist/(maxdist/2)-0.4), 1.0, 0.0])
            cmd.color(col_name, d_name)
            #cmd.color("forest", d_name)
        #elif dist < 10:
            #cmd.color("tv_orange", d_name)
        elif dist > maxdist/2 and dist < maxdist:
            #cmd.set_color(col_name, [0.0, 0.1 + score,1.0 - score])
            cmd.set_color(col_name, [1.0, 1.0 - (dist/(maxdist/2)-1), 0.0])
            cmd.color(col_name, d_name)
        else:
            cmd.color("red", d_name)

        count = count +1

    cmd.set('ray_shadows', 0)
    cmd.set('antialias', 1)
    cmd.ray(1635,1038)
    cmd.png(prot_name + '.png')


if __name__ == "__main__":

    prot_name = sys.argv[1]
    seqfile = sys.argv[2]
    pdbfile = sys.argv[3]
    cfile = sys.argv[4]
    factor = float(sys.argv[5])

    draw(prot_name, seqfile, pdbfile, cfile, factor)
    #draw("1.run_6.S_00000086 1FNAA.fa", "1.run_6.S_00000086.pdb",
    #        "1FNAA.fa.pconsc2.out", 1.0)

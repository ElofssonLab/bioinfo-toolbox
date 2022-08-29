import sys
import operator
import numpy as np
from collections import defaultdict

def printf(format, *args):
        sys.stdout.write(format % args)

def write_pdb_atm_record(atom):
    print("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format("ATOM",atom['atm_no'],atom['atm_name']," ",atom['res_name'],atom['chain'],atom['res_no']," ",atom['x'],atom['y'],atom['z'],atom['occ'],atom['B']))
    return
    
def parse_atm_record(line):

    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    
    return record


def read(pdbfile, chain='', model=1):

    header = ''
    res_lst = []
    atm_lst = []
    tail = ''

    seen_atoms = False
    in_atoms = False
    in_model = False
    curr_resi = 0
    prev_resi = 0
    
    for line in pdbfile:
        if line.startswith('MODEL') and int(line.strip().split()[-1]) == model:
            in_model = True
            header += line
        elif in_model and line.startswith('TER'):
            atm_record = parse_atm_record(atm_lst[-1])
            if chain and not chain == atm_record['chain'] and chain != '*':
                continue
            seen_atoms = True
            #print "seen_atoms model"
            #print len(res_lst)
            #in_atoms = False
            tail += line
        elif in_model and line.startswith('ENDMDL'):
            in_model = False
            tail += line
        elif line.startswith('MODEL') and not int(line.strip().split()[-1]) == model:
            continue
        elif not in_model and line.startswith('TER'):
            atm_record = parse_atm_record(atm_lst[-1])
            if chain and not chain == atm_record['chain']  and chain != '*':
                continue
            seen_atoms = True
            #in_atoms = False
            #print "seen atoms"
            #print len(res_lst)
            continue
        elif not in_model and line.startswith('ENDMDL'):
            continue
        elif not line.startswith('ATOM') and not seen_atoms:
            header += line
        elif not line.startswith('ATOM') and seen_atoms:
            tail += line
        elif in_model or ((not seen_atoms) and (not in_model)):
            atm_record = parse_atm_record(line)
            if chain and not chain == atm_record['chain']  and chain != '*':
                atm_lst = [line]
                continue
            if not in_atoms:
                curr_resi = atm_record['res_no']
                prev_resi = curr_resi
            in_atoms = True
            curr_resi = atm_record['res_no']
            if curr_resi == prev_resi:
                atm_lst.append(line)
            else:
                res_lst.append(atm_lst)
                atm_lst = [line]
            prev_resi = curr_resi
    res_lst.append(atm_lst)
     
    pdbfile.close()
    pdb_lst = [header, res_lst, tail]
    return pdb_lst


def read_chain(pdbfile, chain):

    header = ''
    res_lst = []
    atm_lst = []
    tail = ''

    seen_atoms = False
    curr_resi = 0
    prev_resi = 0
    
    for line in pdbfile:
        # I think we should skip these when we read a chain
        if line.startswith('HETATM'):
            continue
        if not line.startswith('ATOM') and not seen_atoms:
            header += line
        elif not line.startswith('ATOM') and seen_atoms:
            tail += line
        else:
            atm_record = parse_atm_record(line)
            if not atm_record['chain'] == chain:
                continue
            if not seen_atoms:
                curr_resi = atm_record['res_no']
                prev_resi = curr_resi
            seen_atoms = True
            curr_resi = atm_record['res_no']
            if curr_resi == prev_resi:
                atm_lst.append(line)
            else:
                #atm_lst.append(line)
                res_lst.append(atm_lst)
                atm_lst = [line]
            prev_resi = curr_resi
    res_lst.append(atm_lst)
     
    pdbfile.close()
    pdb_lst = [header, res_lst, tail]
    return pdb_lst



def write(pdb_lst, outfile):

    outfile.write(pdb_lst[0])

    for res in pdb_lst[1]:
        for atm in res:
            outfile.write(atm)
            
    outfile.write(pdb_lst[2])
    outfile.close()


def get_coordinates(pdbfile, chain):

    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain  and chain != '*':
            continue

        res_i = atm_record['res_no']

        if res_dict.keys():
            min_res_i = min(res_dict.keys())
        else:
            min_res_i = res_i
        if res_i > 1000 and len(res_dict) < 1000 and min_res_i + len(res_dict) < 1000:
            continue
        atm = [atm_record['x'], atm_record['y'], atm_record['z']]

        """
        line_arr = line.split()
        #print line_arr

        if line_arr[2].startswith('H'):
            continue

        if len(line_arr[2]) > 4:
            if line_arr[3] != chain:
                continue
            try:
                res_i = int(line_arr[4])
            except ValueError as exc:
                continue
            try:
                atm = map(float, line_arr[5:8])
            except ValueError as exc:
                atm = [float('inf'), float('inf'), float('inf')]
        else:
            if line_arr[4] != chain:
                continue
            try:
                res_i = int(line_arr[5])
            except ValueError as exc:
                continue
            try:
                atm = map(float, line_arr[6:9])
            except ValueError as exc:
                atm = [float('inf'), float('inf'), float('inf')]
        """

        res_dict[res_i].append(np.array(atm))
        
    pdbfile.close()
    return sorted(res_dict.items(), key=operator.itemgetter(0))


def get_res_dict(pdbfile, chain):
# This is the orignal one that actually used the resid fiels
    cb_lst = []
    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue

        atm_record = parse_atm_record(line)

        if atm_record['chain'] != ' ' and atm_record['chain'] != chain  and chain != '*':
            continue

        res_i = atm_record['res_no']
        
        if res_dict.keys():
            min_res_i = min(res_dict.keys())
        else:
            min_res_i = res_i

        # I do not really understand whatis tested here, but seems to never be true/AE
        if res_i > 1000 and len(res_dict) < 1000 and min_res_i + len(res_dict) < 1000:
            continue

        # Why ?
        if atm_record['insert'] == 'X':
            res_i = res_i * 0.001

        atm = [float('inf'), float('inf'), float('inf')]

        if atm_record['atm_name'] == 'CA':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm))   
        elif atm_record['atm_name'] == 'CB':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm)) 
    
    return res_dict


def get_bfactor_area(pdbfile, chain):
    # Assumes that the B-factor is used to describe exposed surface area.
    #  Returns surface area and relatice surface are per residue
    
    b_dict = defaultdict(list)
    maxarea={"A":115,"R":225,"D":150,"N":160,"C":135,"E":190,"Q":180,"G":75,"H":195,"I":175,"L":170,"K":200,"M":185,"F":210,"P":145,"S":115,"T":140,"W":255,"Y":230,"V":155}
    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
    #one_to_three = {'R':'ARG', 'H':'HIS', 'K':'LYS', 'D':'ASP', 'E':'GLU', 'S':'SER', 'T':'THR', 'N':'ASN', 'Q':'GLN', 'C':'CYS', 'G':'GLY', 'P':'PRO', 'A':'ALA', 'I':'ILE', 'L':'LEU', 'M':'MET', 'F':'PHE', 'W':'TRP', 'Y':'TYR', 'V':'VAL', 'X': 'UNK'}
    lastres=-9999
    area=0.0
    fracarea=0.0
    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue

        atm_record = parse_atm_record(line)

        if atm_record['chain'] != ' ' and atm_record['chain'] != chain  and chain != '*':
            continue

        res_i = atm_record['res_no']
        
        if b_dict.keys():
            min_res_i = min(b_dict.keys())
        else:
            min_res_i = res_i
        if res_i > 1000 and len(b_dict) < 1000 and min_res_i + len(b_dict) < 1000:
            continue

        if atm_record['insert'] == 'X':
            res_i = res_i * 0.001

        #atm = [float('inf'), float('inf'), float('inf')]
        if res_i != lastres:
            if lastres != -9999 :
                bfactor = [float(area), float(fracarea)]
                #print lastres,bfactor
                b_dict[lastres].append(np.array(bfactor))
            lastres = res_i
            area=0.
            fracarea=0.
        area+=atm_record['B']
        fracarea+=atm_record['B']/maxarea[three_to_one[atm_record['res_name']]]
        
    bfactor = [float(area), float(fracarea)]    
    b_dict[lastres].append(np.array(bfactor))
    return b_dict


def get_area(pdbfile, chain):

    res_dict = get_bfactor_area(pdbfile, chain)

    bfactor_lst = []

    # need to sort to get the sequence correct
    sorted_keys = sorted(res_dict.keys())
    
    for i in sorted_keys:
        bfactor_lst.append(res_dict[i][0])
    pdbfile.close()
    return bfactor_lst


def get_dist_to_surface(pdbfile, chain):

    bf_dict = get_bfactor_area(pdbfile, chain)

    bfactor_lst = []

    # need to sort to get the sequence correct
    sorted_keys = sorted(bf_dict.keys())

    for i in sorted_keys:
        bfactor_lst.append(bf_dict[i][0])

    pdbfile.seek(0,0)

    res_dict = get_res_dict(pdbfile, chain)
    cb_lst = []
    tmp_i = 0

    # need to sort to get the sequence correct
    
    for i in sorted_keys:
        if len(res_dict[i]) > 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][-1])
        elif len(res_dict[i]) == 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][0])
            #print tmp_i,i,res_dict[i][0],res_dict[i][-1]
    pdbfile.close()

    dist_lst = []
    exposed = 0.25
    k=0
    for i in sorted_keys:
        if bfactor_lst[k][1]>exposed:
            dist_lst.append(0.0)
        else:
            dist_lst.append(9999.0)
            l=0
            for j in sorted_keys:
                if k != l:
                    if bfactor_lst[l][1]>exposed:
                        dist_vec=cb_lst[k]-cb_lst[l]
                        dist=np.sqrt(np.sum(dist_vec * dist_vec))
                        if dist < dist_lst[k]:
                            dist_lst[k]=dist

                l+=1
        k+=1
    return dist_lst

def get_ca_coordinates(pdbfile, chain):

    res_dict = get_res_dict(pdbfile, chain)

    ca_lst = []

    # need to sort to get the sequence correct
    sorted_keys = sorted(res_dict.keys())
    
    for i in sorted_keys:
        ca_lst.append(res_dict[i][0])
    pdbfile.close()
    return ca_lst


def get_cb_coordinates(pdbfile, chain):

    res_dict = get_res_dict(pdbfile, chain)

    cb_lst = []
    tmp_i = 0

    # need to sort to get the sequence correct
    sorted_keys = sorted(res_dict.keys())
    
    for i in sorted_keys:
        if len(res_dict[i]) > 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][-1])
        elif len(res_dict[i]) == 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][0])
            #print tmp_i,i,res_dict[i][0],res_dict[i][-1]
    pdbfile.close()
    return cb_lst



def get_atom_seq(pdbfile, chain='', model=1, return_lines=False):

    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
    line_dict = defaultdict(list)
    res_dict = {}
    
    in_model = False

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    #print "Using Chain:", chain
    
    res_name = ''
    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain and chain != '*':
            continue
        res_i = atm_record['res_no']

        if res_dict.keys():
            min_res_i = min(res_dict.keys())
        else:
            min_res_i = res_i
        if res_i > 1000 and len(res_dict) < 1000 and min_res_i + len(res_dict) < 1000:
            continue

        if atm_record['insert'] == 'X':
            res_i = res_i * 0.001
         
        line_dict[res_i].append(line)
        
        if atm_record['atm_name'] != 'CA':
            continue

        if atm_record['res_name'] in three_to_one:
            #res_name = three_to_one[atm_record['res_name']]
            #print res_name
            res_name = three_to_one[atm_record['res_name']]
        #else:
            #res_name = ''
            #continue
        res_dict[res_i] = res_name

    line_lst = [
            l[1]
            for l in sorted(line_dict.items(), key=operator.itemgetter(0))]
    res_lst = sorted(res_dict.items(), key=operator.itemgetter(0))
    atom_seq = ''

    for res in res_lst:
        atom_seq += res[1]

    pdbfile.close()
    if return_lines:
        return atom_seq, ['', line_lst, 'END\n']
    return atom_seq


def get_first_chain(pdbfile):

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        break

    return atm_record['chain']
 
def get_all_chains(pdbfile):
    chains=[]
    lastchain='*'
    for line in open(pdbfile, 'r'):
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != lastchain:
            chains.append(atm_record['chain'])
            lastchain=atm_record['chain']
    return chains
    

def get_acc(pdbfile):

    for line in pdbfile:
        if line.startswith('HEADER'):
            return line[62:66].lower()
    # if no header line in pdb file:
    return ''


if __name__ == '__main__':

    pdbfile = open(sys.argv[1], 'r')
    if len(sys.argv) == 3:
        chain = sys.argv[2]
    else:
        chain = get_first_chain(pdbfile)
    print(get_atom_seq(pdbfile, chain))
    pdbfile.close()
    #pdbfile = open(sys.argv[1], 'r')
    #print get_coordinates(pdbfile)
    #pdbfile.close()

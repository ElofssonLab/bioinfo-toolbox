import sys
from collections import defaultdict


# read rosetta vall fragment database and parse it into:
# [header_lines, {pdb_chain_id : [fragments]}]
# with:
#  - header_lines: [lines starting with '#']
#  - pdb_chain_id: '4letter-PDB' + 'chain', eg: 1xxxA
#  - fragments: list of lines starting with pdb_chain_id
def read(vallfile):

    header = []
    id_list = []
    frag_dict = defaultdict(list)

    for line in vallfile:
        if line.startswith('#'):
            header.append(line)
        else:
            id = line.split()[0][:4]
            frag_dict[id].append(line)
            if not id in id_list:
                id_list.append(id)

    vallfile.close()
    result = [header, id_list, frag_dict]
    return result


# write [header_lines, {pdb_chain_id : [fragments]}]
# into rosetta vall fragment database format
def write(vallfile_list, outfile):

    header = vallfile_list[0]
    id_list = vallfile_list[1]
    frag_dict = vallfile_list[2]
    
    for line in header:
        outfile.write(line)

    for id in id_list:
        if id in frag_dict:
            for line in frag_dict[id]:
                outfile.write(line)

    outfile.close()


if __name__ == '__main__':
    
    vallfile_list = read(open(sys.argv[1], 'r'))
    write(vallfile_list, open('%s.out' % sys.argv[1], 'w'))





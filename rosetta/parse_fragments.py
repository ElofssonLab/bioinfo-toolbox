import sys
from collections import defaultdict


# read rosetta fragment file and parse it into:
# [(position_header, fragment_dictionary)]
def read(fragfile):

    fragfile_list = []
    frag_dict = defaultdict(list)

    key = ''
    last_was_header = False

    for line in fragfile:

        if line.strip() == '':
            continue

        if line.startswith('position:') and key != '':
            fragfile_list.append((key, frag_dict))
            frag_dict = defaultdict(list)
            key = line

        elif line.startswith('position:') and key == '':
            key = line

        elif not line.startswith('position:'):
            line_arr = line.split()
            id = line.split()[0]
            frag_dict[id].append(line)
            last_was_header = False

    fragfile_list.append((key, frag_dict))
    fragfile.close()
    return fragfile_list


# write [(position_header, fragment_dictionary)]
# into rosetta fragment file format
def write(fragfile_list, outfile, frag_len):

    for (key, frag_dict) in fragfile_list:
        outfile.write('%s\n' % key)
        for id, frag_lines in frag_dict.iteritems():
            count = 1
            for line in frag_lines:
                outfile.write(line)
                if count % frag_len == 0:
                    outfile.write('\n')
                count += 1

    outfile.close()


# return length of the fragments in current fragment file
def get_frag_len(fragfile):

    frag_len = 0
    last_was_header = False

    for line in fragfile:

        if line.startswith('position:'):
            last_was_header = True
        elif (not line.startswith('position:')) and line.strip() != '':
            last_was_header = False
            frag_len += 1
        elif line.strip() == '' and not last_was_header:
            break

    fragfile.close()
    return frag_len


if __name__ == '__main__':
    
    fragfile_list = read(open(sys.argv[1], 'r'))

    print len(fragfile_list)
    print 'fragment length: %d' % get_frag_len(open(sys.argv[1], 'r'))

    write(fragfile_list, open('%s.out' % sys.argv[1], 'w'), 3)





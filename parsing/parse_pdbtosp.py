import sys

def parse_line(line):
    
    pdb_id = line[:6].strip()
    method = line[6:16].strip()
    resolu = line[16:28].strip()
    sp_id = line[28:40].strip()

    line_arr = [pdb_id, method, resolu, sp_id]

    return line_arr


def parse(infile):

    result_dict = {}
    in_map = False

    for line in infile:

        if line.startswith("____"):
            in_map = True
            continue
        if not line.strip():
            in_map = False
        if in_map:
            line_arr = parse_line(line.strip())
            result_dict[line_arr[0]] = (line_arr[1], line_arr[2], line_arr[3])
            continue
    return result_dict


if __name__ == "__main__":
    
    mapping = parse(open(sys.argv[1], 'r'))
    print len(mapping)

    print mapping["1A3A"]


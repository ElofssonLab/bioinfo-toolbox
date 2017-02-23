import sys
from collections import defaultdict

def read(f):
    """ read SEG output and store it in:
        {parameters: {id: sequence}}
        where lower-case letters in "sequence" denote low-complexity
        under given parameters
    """
    result_dict = defaultdict(dict)
    key = ''
    seq = ''
    for l in f:
        if 'low complexity' in l:
            params = l.split(':')[-1].strip()
            continue
        if l.startswith('>'):
            # not first entry
            if key:
                result_dict[params][key] = seq
                seq = ''
            key = l.split('|')[1]
            continue
        # skip empty lines
        if not l.strip():
            continue
        # lower letters to the left
        if l.strip()[0].islower():
            seq += l.strip().split()[0]
        # capitals to the right
        else:
            seq += l.strip().split()[-1]
    return result_dict


if __name__ == "__main__":

    print read(open(sys.argv[1]))
                

        

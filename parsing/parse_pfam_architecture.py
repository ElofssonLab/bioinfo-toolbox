import sys
from collections import defaultdict

def get_arch_dict(infile):
    
    result = defaultdict(list)
    for l in infile:
        l = l.strip()
        l_arr = l.split('\t')
        arch_name = l_arr[1]
        dom_arr = l_arr[-1].split(' ')
        for dom in dom_arr:
            if not arch_name in result[dom]:
                result[dom].append(arch_name)
    return result


if __name__ == '__main__':

    print len(get_arch_dict(open(sys.argv[1])))

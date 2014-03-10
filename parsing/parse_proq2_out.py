import sys

def parse(f):
    result = []
    f.readline() # skip line header
    for l in f:
        l = l.strip()
        l_lst = l.split()
        map(float, l_lst[1:-1])
        result.append(l_lst[1:])
    # sort by overall score
    result.sort(key=lambda itm: itm[0], reverse=True)
    return result


if __name__ == "__main__":
    result = parse(open(sys.argv[1], 'r'))
    for sc in result:
        print ' '.join(sc)

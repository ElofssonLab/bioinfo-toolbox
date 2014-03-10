import sys

def parse(f, prefix=""):
    result = []
    f.readline() # skip line header
    for l in f:
        l = l.strip()
        l_lst = l.split()
        map(float, l_lst[1:-1])
        new_tag = prefix + l_lst[-1]
        l_lst[-1] = new_tag
        result.append(l_lst[1:])
    # sort by overall score
    result.sort(key=lambda itm: itm[0], reverse=True)
    return result


if __name__ == "__main__":
    if len(sys.argv) == 3:
        prefix = sys.argv[2]
    else:
        prefix = ""
    result = parse(open(sys.argv[1], 'r'), prefix)
    for sc in result:
        print ' '.join(sc)

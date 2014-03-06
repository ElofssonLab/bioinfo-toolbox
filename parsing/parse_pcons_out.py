import sys

def parse(f):
    result = []
    in_flag = False
    for l in f:
        if l.startswith("QMODE"):
            in_flag = True
            continue
        if l.startswith("END"):
            in_flag = False
            continue
        if in_flag:
            l_lst = l.split()
            result.append((l_lst[0], float(l_lst[1])))
    result.sort(key=lambda itm: itm[1], reverse=True)
    return result


if __name__ == "__main__":
    
    # Parse single file and print full list 
    if len(sys.argv) == 2:
        f = open(sys.argv[1],'r')
        scores = parse(f)
        f.close()
        for sc in scores:
            print ' '.join(map(str,sc))
    # Parse multiple files and print top-ranked
    elif len(sys.argv) > 2:
        for fname in sys.argv[1:]:
            f = open(fname,'r')
            scores = parse(f)
            f.close()
            print ' '.join(map(str,scores[0]))


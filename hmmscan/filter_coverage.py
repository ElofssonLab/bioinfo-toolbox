import sys


def filter_coverage(lines, th=0.9):
    result = []
    for l in lines:
        if l.startswith('#'):
            result.append(l)
            continue
        l_arr = l.strip().split()
        length = float(l_arr[2])
        start = float(l_arr[15])
        end = float(l_arr[16])
        if (end-start) / length > th:
            result.append(l)
    return result


if __name__ == '__main__':

    fname = sys.argv[1]
    if len(sys.argv) == 3:
        th = float(sys.argv[2])
    else:
        th = 0.9

    with open(fname) as f:
        lines = f.readlines()
    lines_filtered = filter_coverage(lines, th=th)
    for l in lines_filtered:
        sys.stdout.write(l)


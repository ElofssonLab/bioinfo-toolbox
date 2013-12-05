import sys
from operator import itemgetter

def parse(infile):
    result = []
    in_summary = False
    for line in infile:
        if line.startswith(">> Summary "):
            in_summary = True
            continue
        if in_summary:
            if line.startswith("Filename"):
                continue
            elif line.startswith("---------"):
                continue
            elif not line.strip():
                break
            else:
                line_arr = line.split()
                # entry = (filename, molpdf, DOPE, GA341)
                result.append((line_arr[0], float(line_arr[1]), float(line_arr[2]), float(line_arr[3])))

    infile.close()
    return sorted(result, key=itemgetter(2))


if __name__ == "__main__":
    print parse(open(sys.argv[1], 'r'))





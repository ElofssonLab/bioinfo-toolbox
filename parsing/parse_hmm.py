import sys

def read_hmm(afile):

    hmm_list = []
    null_model = []
    in_hmm = False

    for line in afile:
        line = line.strip()
        if line.startswith('NULL'):
            null_model = line.split()[1:]
            print null_model

        if line.startswith('HMM'):
            in_hmm = True
            continue

        if in_hmm:
            line_arr = line.split()
            if len(line_arr) != 23:
                continue
            aa = line_arr[0]
            freq = line_arr[2:]
            hmm_list.append((aa, freq))
            
    return hmm_list


if __name__ == "__main__":

    afile = open(sys.argv[1], 'r')
    hmm_list = read_hmm(afile)
    afile.close()

    #print hmm_list

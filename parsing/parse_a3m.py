import sys


def reformat(seq_0, seq_1):
    result_0 = ''
    result_1 = ''
    assert len(seq_0) <= len(seq_1)

    i = 0
    j = 0
    while i < len(seq_0):
        aa_0 = seq_0[i]
        aa_1 = seq_1[j]
        if aa_1.islower():
            result_0 += '-'
            result_1 += aa_1.upper()
            j += 1
        else:
            result_0 += aa_0
            result_1 += aa_1
            i += 1
            j += 1

    return result_0, result_1



def get_pairwise(ali_file, acc_a):
    with open(ali_file) as f:
        seq_0 = ''
        seq_0_tmp = ''
        seq_1 = ''
        in_seq_1 = False
        for l in f:
            if l.startswith('#'):
                continue
            if not seq_0 and not l.startswith('>'):
                seq_0_tmp += l.strip()
            elif not seq_0 and l.startswith('>'):
                seq_0 = seq_0_tmp

            if l.startswith('>'):
                in_seq_1 = acc_a in l
                continue
            if in_seq_1:
                seq_1 += l.strip()
        seq_0, seq_1 = reformat(seq_0, seq_1)
        return seq_0, seq_1
    

if __name__ == '__main__':

    ali_file = sys.argv[1]
    acc = sys.argv[2]
    print get_pairwise(ali_file, acc)

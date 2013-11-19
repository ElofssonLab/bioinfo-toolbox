import sys


def get_scores(bfile):

    result_dict = {}
    ID = ''
    eval = -1.0
    score = -1.0
    ident = -1.0
    posit = -1.0
    len = -1
    
    for line in bfile:
       
        line = line.strip()

        # found new entry: fininsh old and create new
        if line.startswith('>'):
            if ID:
                #score_tup = (eval, ident, posit, score, len)
                score_tup = (eval, ident, score, len)
                result_dict[ID] = score_tup
            ID = line.split(' ')[0]
            ID = ID[1:]

        elif line.startswith('Score = '):
            line_arr = line.split(',')
            score = float(line_arr[0].strip().split(' ')[-3])
            eval = line_arr[1].strip().split(' ')[2]
        elif line.startswith('Identities = '):
            line_arr = line.split(',')
            ident = float(line_arr[0].strip().split(' ')[3][1:-2]) / 100.0
            posit = float(line_arr[1].strip().split(' ')[3][1:-2])
            len = int(line_arr[0].split('/')[1].split(' ')[0])

    # add last entry
    score_tup = (eval, ident, score, len)
    result_dict[ID] = score_tup

    bfile.close()
    return result_dict


if __name__ == '__main__':

    bfile = open(sys.argv[1], 'r')
    print get_scores(bfile)

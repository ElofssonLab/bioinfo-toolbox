import sys


def get_scores(f):
    
    result_dict = {}
    ID = ''
    eval = -1.0
    score = -1.0
    ident = -1.0
    len = -1
    
    in_hit = False

    for line in f:
       
        line = line.strip()

        # look for hit list start
        if line.startswith('No Hit'):
            in_hit = True
            continue

        # end of hit list        
        if not line and in_hit:
            in_hit = False
            ID = ''
            continue

        # parse hit list
        if in_hit:
            # get rid of 'HMM' column
            line_arr_tmp = line.split('(')
            line_tmp = ''.join(line_arr_tmp[:-1])
            line_arr = line_tmp.split()
            #line_arr = line.split()
            ID = line_arr[1]
            eval = line_arr[-7]
            score = float(line_arr[-5])
            len = int(line_arr[-3])
            # the only thing missing is ident (default: -1)
            if not ID in result_dict:
                result_dict[ID] = (eval, ident, score, len)
            continue

        # parse verbose result list to get the identities
        if line.startswith('>'):
            if ID and result_dict[ID][1] == -1:
                score_lis = list(result_dict[ID])
                score_lis[1] = ident
                result_dict[ID] = tuple(score_lis)
                #print score_lis
            ID = line.split(' ')[0]
            ID = ID[1:]
            #print ID
        elif line.startswith('Probab='):
            line_arr = line.split()
            ident = float(line_arr[4].strip().split('=')[1][:-1]) / 100.0

    # update last entry
    score_lis = list(result_dict[ID])
    score_lis[1] = ident
    result_dict[ID] = tuple(score_lis)

    f.close()
    return result_dict


if __name__ == '__main__':

    bfile = open(sys.argv[1], 'r')
    print get_scores(bfile)

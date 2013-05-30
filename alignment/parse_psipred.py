def horizontal(afile):

    result = ''
    for line in afile:
        line_arr = line.strip().split(' ')
        if line_arr[0] == 'Pred:':
            result += line_arr[1]

    return result


def horizontal_seq(afile):

    result = ''
    for line in afile:
        line_arr = line.strip().split(' ')
        if line_arr[0] == 'AA:':
            result += line_arr[1]

    return result

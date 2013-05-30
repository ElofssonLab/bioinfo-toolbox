def get_resolutions():

    resolu_path = '/bubo/home/h9/mircomic/glob/databases/PDB/current_release/resolu.idx'
    resfile = open(resolu_path,'r')
    in_idx = False

    result_dict = {}

    for line in resfile:

        if in_idx:
            line_arr = line.strip().split('\t')
            id = line_arr[0].lower()
            resolu = line_arr[-1]
            if resolu != ';':
                result_dict[id] = float(resolu)

        if line.startswith('--'):
            in_idx = True

    return result_dict


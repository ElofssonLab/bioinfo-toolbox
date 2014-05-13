import sys

def parse(f):
    """ Parses MolProbity oneline output:
    #pdbFileName:x-H_type:chains:residues:nucacids:resolution:rvalue:rfree:clashscore:clashscoreB<40:minresol:maxresol:n_samples:pct_rank:pct_rank40:cbeta>0.25:numCbeta:rota<1%:numRota:ramaOutlier:ramaAllowed:ramaFavored:numRama:numbadbonds:numbonds:pct_badbonds:pct_resbadbonds:numbadangles:numangles:pct_badangles:pct_resbadangles:MolProbityScore:Mol_pct_rank
    """
    
    result = {}
    avg = []
    n = 0

    for l in f:
        # skip header/comments
        if l.startswith('#'):
            continue
        n += 1
        l_lst = l.strip().split(':')
        res_line = []
        for i in xrange(len(l_lst)):
            try:
                res_line.append(float(l_lst[i]))
            except:
                res_line.append(l_lst[i])
        if not avg:
            avg = res_line
            avg[0] = 'Mean'
        else:
            for i in xrange(1, len(res_line)):
                if isinstance(res_line[i], float):
                    if avg[i]:
                        avg[i] += res_line[i]
                    else:
                        avg[i] = res_line[i]
        result[res_line[0]] = res_line[1:]
    for i in xrange(len(avg)):
        if isinstance(avg[i], float):
            avg[i] /= n
    result[avg[0]] = avg
    return result      



if __name__ == "__main__":
    
    result =  parse(open(sys.argv[1], 'r'))

    print result['Mean']


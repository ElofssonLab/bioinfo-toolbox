import sys
import os

IGNORE_LST = ['PFRMAT','TARGET','AUTHOR','REMARK','METHOD','MODEL','QMODE','END']

def parse_name(model):
    model_arr = model.split('.')
    acc = '.'.join(model_arr[:2])
    no = model_arr[3].strip('fa_')
    method = '.'.join(model_arr[4:6])
    stage = model_arr[6]
    return "%s %s %s %s" % (acc, no, method, stage)


def parse(fname):
    with open(fname) as f:
        model = ""
        score_lst = []
        for l in f:
            l_arr = l.strip().split()
            if not l_arr:
                continue
            try:
                first = l_arr[0]
            except IndexError:
                raise IndexError("%s %s" % (fname, l))
            if first in IGNORE_LST or not '.pdb' in first:
                continue
            model = parse_name(os.path.basename(first))
            score_lst = l_arr[2:]
            for i, score in enumerate(score_lst):
                print "%s %d %s" % (model, i+1, score)


if __name__=="__main__":

    fname = sys.argv[1]
    parse(fname)

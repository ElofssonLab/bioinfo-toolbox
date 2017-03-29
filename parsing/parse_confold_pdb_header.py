import sys
import os

def parse_name(model):
    model_arr = model.split('.')
    acc = '.'.join(model_arr[:2])
    no = model_arr[-5].strip('fa_')
    method = '.'.join(model_arr[-4:-2])
    stage = model_arr[-2]
    return "%s %s %s %s" % (acc, no, method, stage)


def parse(flist):
    for fname in flist:
        with open(fname) as f:
            model = parse_name(os.path.basename(fname))
            for l in f:
                if not l.startswith('REMARK'):
                    break
                l_arr = l.split()
                if l_arr[-1].startswith('FILENAME'):
                    orig_fname = l_arr[-1].split('=')[-1].strip('"')
                    model_raw = orig_fname.split('.')[-2].strip('fa_')
                if l.startswith('REMARK overall'):
                    total = l_arr[-1]
                if l.startswith('REMARK bon'):
                    bond = l_arr[-1]
                if l.startswith('REMARK ang'):
                    angle = l_arr[-1]
                if l.startswith('REMARK imp'):
                    imp = l_arr[-1]
                if l.startswith('REMARK vdw'):
                    vdw = l_arr[-1]
                if l.startswith('REMARK noe'):
                    noe = l_arr[-1]
            print "%s %s %s %s %s %s %s %s" % (model, model_raw, total, bond, angle, imp, vdw, noe)


if __name__=="__main__":

    flist = sys.argv[1:]
    parse(flist)

import sys

IGNORE_LST = ['PFRMAT','TARGET','AUTHOR','REMARK','METHOD','MODEL','QMODE','END']

def reformat(fname):
    with open(fname) as f:
        for l in f:
            l_arr = l.strip().split()
            first = l_arr[0]
            if first in IGNORE_LST or not '.pdb' in first:
                continue
            first_arr = first.split('.')
            if first_arr[-4].startswith('confold'):
                setting = '.'.join(first_arr[-4:-2])
                model = first_arr[-5].split('_')[-1]
            else:
                setting = first_arr[-3]
                model = first_arr[-4].split('_')[-1]
            stage = first_arr[-2] 

            name = '.'.join(first_arr[:2])

            print name, model, setting, stage, l_arr[1]



if __name__ == '__main__':
    reformat(sys.argv[1])

import sys
from math import *
from collections import defaultdict

import matplotlib.pyplot as plt


def s2d(s):
    d0=sqrt(5)
    d=100 # for CASP we cap the distance at 100 angstroms
    if s>0.0004: # this is the S score for 100 angstroms
        if s>=1:
            d=0
        else:
            d=sqrt(1/s-1)*d0
    return d 


def parse_qa(f):
    qa_dict = defaultdict(list)
    ignore_lst = ['PFRMAT', 'TARGET', 'AUTHOR', 'REMARK', 'MODEL', 'QMODE', 'END']
    with open(f) as qa_f:
        acc = ''
        method = ''
        for l in qa_f:
            l_arr = l.strip().split()
            # replace unpredicted positions
            l_arr = [100 if x=='X' else x for x in l_arr]
            if l.startswith('METHOD'):
                if not method:
                    method = l_arr[1]
                    if '*' in method:
                        method = 'Pcons'
                    continue
                else:
                    continue
            if l_arr[0] in ignore_lst:
                continue
            try:
                float(l_arr[0])
                #if method == 'Pcons':
                #    qa_dict[acc] += map(s2d, map(float, l_arr))
                #else:
                qa_dict[acc] += map(float, l_arr)
            except ValueError:
                acc = l_arr[0]
                #if method == 'Pcons':
                #    qa_dict[acc] = map(s2d, map(float, l_arr[2:]))
                #else:
                qa_dict[acc] = map(float, l_arr[2:])
    return method, qa_dict


def plot_qa(f_lst, outf_prefix=''):

    data = defaultdict(dict)
    for f in f_lst:
        method, qa_dict = parse_qa(f)
        for acc, qa_lst in qa_dict.iteritems():
            data[acc][method] = qa_lst

    for acc, method_dict in data.iteritems():
        fig = plt.figure(figsize=(10,4))
        ax = plt.axes()
        plot_lst = []
        label_lst = []
        for method, qa_lst in method_dict.iteritems():
            x = xrange(len(qa_lst))
            plt.plot(x, qa_lst, label=method)
            #plot_lst.append(plt.plot(x, qa_lst))
            #label_lst.append(method)
        #plt.legend(plot_lst, label_lst)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
        #plt.show()
        plt.ylim(0,30)
        outf = '%s%s.qa.png' % (outf_prefix, acc)
        plt.savefig(outf)
        plt.close()

            
            


if __name__ == '__main__':

    outf_prefix = sys.argv[1]
    f_lst = sys.argv[2:]
     
    plot_qa(f_lst, outf_prefix=outf_prefix)


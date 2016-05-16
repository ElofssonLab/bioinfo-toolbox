import sys
from collections import defaultdict

import matplotlib.pyplot as plt


def parse_qa(f):
    qa_dict = defaultdict(list)
    ignore_lst = ['PFRMAT', 'TARGET', 'AUTHOR', 'REMARK', 'MODEL', 'QMODE', 'END']
    with open(f) as qa_f:
        acc = ''
        for l in qa_f:
            l_arr = l.strip().split()
            # replace unpredicted positions
            l_arr = [100 if x=='X' else x for x in l_arr]
            if l.startswith('METHOD'):
                method = l_arr[1]
                continue
            if l_arr[0] in ignore_lst:
                continue
            try:
                float(l_arr[0])
                qa_dict[acc] += map(float, l_arr)
            except ValueError:
                acc = l_arr[0]
                qa_dict[acc] = map(float, l_arr[2:])
    return method, qa_dict


def plot_qa(f_lst, outf_prefix=''):

    data = defaultdict(dict)
    for f in f_lst:
        method, qa_dict = parse_qa(f)
        for acc, qa_lst in qa_dict.iteritems():
            data[acc][method] = qa_lst

    for acc, method_dict in data.iteritems():
        fig = plt.figure(figsize=(10,2))
        ax = plt.axes()
        plot_lst = []
        label_lst = []
        for method, qa_lst in method_dict.iteritems():
            x = xrange(len(qa_lst))
            plt.plot(x, qa_lst, label=method)
            #plot_lst.append(plt.plot(x, qa_lst))
            #label_lst.append(method)
        #plt.legend(plot_lst, label_lst)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)
        #plt.show()
        plt.ylim(0,25)
        outf = '%s%s.qa.png' % (outf_prefix, acc)
        plt.savefig(outf)

            
            


if __name__ == '__main__':

    outf_prefix = sys.argv[1]
    f_lst = sys.argv[2:]
     
    plot_qa(f_lst, outf_prefix=outf_prefix)


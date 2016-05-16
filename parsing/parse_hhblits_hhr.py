import string, copy
import sys
from collections import defaultdict


# returns 'True' for two overlapping intervals, 'False' otherwise
def overlapcheck(pos_x, pos_y):
    
    pos_x_vec = pos_x.split('-')
    pos_y_vec = pos_y.split('-')
    pos_x0 = int(pos_x_vec[0])
    pos_x1 = int(pos_x_vec[1])
    pos_y0 = int(pos_y_vec[0])
    pos_y1 = int(pos_y_vec[1])
    overlap =  max(0, min(pos_x1, pos_y1) - max(pos_x0, pos_y0))

    return overlap != 0


## parses '*.hhr' output from hhblits
def main(afile):

    eval_dict = {}
    pos_dict = {}
    sorted_out_dict = {}
    hmm_count = 0

    in_reslist = False

    for aline in afile:
        aline = aline.strip()

        # have we seen the header yet?
        if in_reslist:
            if not aline == '':
                alist = aline.split()
                key = alist[1] 
                hmm_count += 1
                
                # in some cases the last two elements are not separated by space
                if not '-' in alist[-1]:
                    pos = alist[-2]
                    eval = float(alist[-8])
                else:
                    pos = alist[-1].split('(')[0]
                    eval = float(alist[-7])
                
                # now we have proper pos/e-values
                if eval_dict.has_key(key):

                    # insert only if not overlapping or 
                    # if overlapping and higher evalue
                    insert_flag = True
                    for i in range(len(eval_dict[key])):
                        is_overlap = overlapcheck(pos,pos_dict[key][i])
                        if is_overlap:
                            if eval < eval_dict[key][i]:
                                eval_dict[key].remove(i)
                                pos_dict[key].remove(i)
                            else:
                                insert_flag = False
                    if insert_flag:
                        eval_dict[key].append(eval)
                        pos_dict[key].append(pos)
                    else:
                        if sorted_out_dict.has_key(key):
                            sorted_out_dict[key].append(eval_dict[key][i])
                        else:
                            sorted_out_dict[key] = [eval_dict[key][i]]
                else:
                    eval_dict[key] = [eval]
                    pos_dict[key] = [pos]
            
            # result summary ends with empty line
            else:
                in_reslist = False

        # recognize header of the result summary
        if aline.startswith('No Hit'):
            in_reslist = True

    return pos_dict, eval_dict, sorted_out_dict, hmm_count


def parse_hit_table(afile, th_e=0.0001, th_cov=0.9, th_prob=99.):
    
    hit_dict={}
    in_reslist = False
    cols_qry = 0

    for l in afile:
        l = l.strip()

        if l.startswith('Match_columns'):
            cols_qry = float(l.split()[1])
        
        # this is where we want to be
        if in_reslist:
            if l:
                l_arr = l.split()
                key = l_arr[1] 
            
                # in some cases the last two elements are not separated by space
                if not '-' in l_arr[-1]:
                    prob = float(l_arr[-9])
                    e_val = float(l_arr[-8])
                    p_val = float(l_arr[-7])
                    score = float(l_arr[-6])
                    ss = float(l_arr[-5])
                    cols = float(l_arr[-4])
                    pos_qry = map(int, l_arr[-3].split('-'))
                    pos_tmpl = map(int, l_arr[-2].split('-'))
                else:
                    prob = float(l_arr[-8])
                    e_val = float(l_arr[-7])
                    p_val = float(l_arr[-6])
                    score = float(l_arr[-5])
                    ss = float(l_arr[-4])
                    cols = float(l_arr[-3])
                    pos_qry = map(int, l_arr[-2].split('-'))
                    pos_tmpl = map(int, l_arr[-1].split('(')[0].split('-'))

                cov = cols/cols_qry

                if e_val <= th_e and prob >= th_prob and cov >= th_cov:
                    hit_dict[key] = [prob, e_val, p_val, score, ss, cols, pos_qry, pos_tmpl]
                else:
                    continue
            else:
                in_reslist = False
                break

        # recognize header of the result summary
        if l.startswith('No Hit'):
            in_reslist = True

    return hit_dict


def parse_alignments(pdb_aligment_path, No=-1):

    with open(pdb_aligment_path) as afile:
        ali_dict=defaultdict(list)
        t_id = ""
        t_seq = ""
        q_seq = ""

        in_ali = False
        in_no = False
        for l in afile:
            l = l.strip()
            if l.startswith("No"):
                if No == -1:
                    in_no = True
                else:
                    in_no = str(No) in l
                continue
            if in_no and l.startswith(">"):
                if t_seq and q_seq:
                    ali_dict[t_id].append((t_seq, q_seq))
                t_id = l.strip('>').split()[0]
                q_seq = ""
                t_seq = ""
                if not in_ali:
                    in_ali = True
                    continue
                else:
                    continue
            if in_no and in_ali:
                if l.startswith("Q") and not "Consensus" in l:
                    q_seq += l.split()[3]
                if l.startswith("T") and t_id in l:
                    t_seq += l.split()[3]
    return ali_dict


if __name__ == "__main__":
    #afile = open(sys.argv[1], 'r')
    #pos_dict, eval_dict, sorted_out_dict, hmm_count = main(afile)
    #hit_dict = parse_hit_table(afile)
    #afile.close()

    print parse_alignments(sys.argv[1])

    """
    #print pos_dict
    print hmm_count
    print len(eval_dict)
    print len(sorted_out_dict)
    print float(len(sorted_out_dict))/float(hmm_count)
    """
    for key, hit in hit_dict.iteritems():
        print key, hit
    print "Found %s hits that pass thresholds." % len(hit_dict)

import sys
import string
import csv

#import parse_dom_assign_hmmer3_repeat
import parse_hmmer3_domtblout
#import dict_pfam_pdb as pfam_pdb
#import dict_pfam_pdb_clan as pfam_pdb
#import dict_pfam_numseq as pfam_numseq
#import dict_pfam_numseq_clan as pfam_numseq

def create_dom_dict(domfile, pfam_clan, clan_name, pfameval, extraeval, pdb_res_dict):

    ## parse hmmscan result file into two dictionaries
    #scan_result = parse_dom_assign_hmmer3_repeat.main(domfile, pfam_clan, clan_name, pfameval, extraeval)
    repeat_flag = True
    scan_result = parse_hmmer3_domtblout.main(domfile, pfameval, extraeval, repeat_flag, pfam_clan, clan_name)
    domkey_name_dict = scan_result[0]
    dom_prot_dict = scan_result[1]
    prot_dom_dict = scan_result[2]
    domfile.close()

    #prot_dom_dict = scan_result[0]
    #domkey_name_dict = scan_result[1]
    
    ## get domain statistics
    dom_dict = {}  # [len, avg E-val, avg Bit-score, max Bit-score, avg #rep, min #rep, max #rep, #seq, best PDB, longest PDB, name]
    
    for prot_key, dom_list  in prot_dom_dict.iteritems():
        dom_counts_dict = {}
        dom_pos_dict = {}
        dom_pdb = prot_key # since hmmscan was run on PDB sequences
                
        if len(dom_list) > 0:

            for dom in dom_list:
                dom_pos = dom[0]
                dom_key = dom[2]
                if dom_pos_dict.has_key(dom_key):
                    dom_pos_dict[dom_key].append(dom_pos)
                    dom_pos_dict[dom_key].sort()
                else:
                    dom_pos_dict[dom_key] = [dom_pos]
                if dom_counts_dict.has_key(dom_key):
                    dom_counts_dict[dom_key] += 1
                else:
                    dom_counts_dict[dom_key] = 1


            for dom in dom_list:
                if pdb_res_dict.has_key(dom_pdb[0:4]):
                    dom_pos = dom[0]
                    #dom_eval = dom[1] # not used here
                    dom_key = dom[2]
                    #dom_len = dom[3]
                    dom_len = dom_pos[1] - dom_pos[0]
                    dom_score = dom[3]
                    dom_count = dom_counts_dict[dom_key]
                    dom_res = pdb_res_dict[dom_pdb[0:4]]
                    
                    dom_clan = 'None'
                    if pfam_clan.has_key(string.split(dom_key,'.')[0]):
                        dom_clan = pfam_clan[string.split(dom_key,'.')[0]]
                    dom_name = domkey_name_dict[dom_key]

                    if dom_dict.has_key(dom_key):
                        dom_dict[dom_key][0] = max(dom_dict[dom_key][0], dom_len)
                        #dom_dict[dom_key][2] = (dom_dict[dom_key][2] + dom_score) / 2.0

                        # select best structure
                        old_res = pdb_res_dict[dom_dict[dom_key][6][0:4]]
                        old_score = dom_dict[dom_key][2]
                        if dom_score > old_score:
                            dom_dict[dom_key][2] = dom_score 
                        if dom_res != -1 and old_res != -1:
                            #if (dom_score / dom_res) < (old_score / old_res):
                            if dom_res < 2.8 and dom_score > old_score:
                                dom_dict[dom_key][1] = dom_score 
                                dom_dict[dom_key][6] = dom_pdb
                                dom_dict[dom_key][7] = dom_res
                                dom_dict[dom_key][8] = '%d-%d' % (dom_pos[0], dom_pos[1])
                        elif dom_res != -1:
                            dom_dict[dom_key][1] = dom_score 
                            dom_dict[dom_key][6] = dom_pdb
                            dom_dict[dom_key][7] = dom_res
                            dom_dict[dom_key][8] = '%d-%d' % (dom_pos[0], dom_pos[1])
                        elif dom_res == -1 and old_res == -1:
                            if dom_score > old_score:
                                dom_dict[dom_key][1] = dom_score 
                                dom_dict[dom_key][6] = dom_pdb
                                dom_dict[dom_key][7] = dom_res
                                dom_dict[dom_key][8] = '%d-%d' % (dom_pos[0], dom_pos[1])

                        # select structure with most repeats
                        old_max_count = dom_dict[dom_key][5]
                        #old_min_pos = int(dom_dict[dom_key][11].split('-')[0])
                        #old_max_pos = int(dom_dict[dom_key][11].split('-')[1])
                        dom_dict[dom_key][3] = (dom_dict[dom_key][3] + dom_count) / 2
                        dom_dict[dom_key][4] = min(dom_dict[dom_key][4], dom_count)
                        dom_dict[dom_key][5] = max(old_max_count, dom_count)
                        if dom_count > old_max_count:
                            dom_dict[dom_key][9] = dom_pdb
                            dom_dict[dom_key][10] = dom_res
                            dom_dict[dom_key][11] = dom_pos_dict[dom_key]
                            #dom_dict[dom_key][11] = '%d-%d' % (min(old_min_pos, dom_pos[0]),
                            #                                   max(old_max_pos, dom_pos[1]))
                        elif dom_count == old_max_count:
                            if dom_res < old_res:
                                dom_dict[dom_key][9] = dom_pdb
                                dom_dict[dom_key][10] = dom_res
                                dom_dict[dom_key][11] = dom_pos_dict[dom_key]
                                #dom_dict[dom_key][11] = '%d-%d' % (min(old_min_pos, dom_pos[0]),
                                #                                   max(old_max_pos, dom_pos[1]))


                    else:
                        dom_dict[dom_key] = [dom_len,
                                            dom_score, # score in best PDB
                                            dom_score, # max score
                                            dom_count, # avg
                                            dom_count, # min
                                            dom_count, # max
                                            dom_pdb, # best PDB
                                            pdb_res_dict[dom_pdb[0:4]], # res in best PDB
                                            '%d-%d' % (dom_pos[0], dom_pos[1]),
                                            dom_pdb, # longest PDB
                                            pdb_res_dict[dom_pdb[0:4]], # res in longest PDB
                                            #'%d-%d' % (dom_pos[0], dom_pos[1]),
                                            dom_pos_dict[dom_key],
                                            dom_clan,
                                            dom_name]
            
    return dom_dict


def write_dom_dict(dom_dict, outfile):
    #sorted(pfam_stat_dict.items(), key = lambda x : x[1][1], reverse = True)[0:10]

    writer = csv.writer(outfile, delimiter='\t')
    colnames = ['PfamID', 'length', 'score_in_best_PDB', 'max_score', 'avg_num_repeats', 'min_num_repeats', 'max_num_repeats', 'best_PDB', 'res_in_best_PDB', 'pos_in_best_PDB', 'longest_PDB', 'res_in_longest_PDB', 'pos_in_longest_PDB', 'ClanID', 'Name']
    writer.writerow(colnames)

    for dom_key, dom_statlist in dom_dict.iteritems():
        line = [dom_key] + dom_statlist
        #print [dom_key] + dom_statlist
        #print line
        writer.writerow(line)


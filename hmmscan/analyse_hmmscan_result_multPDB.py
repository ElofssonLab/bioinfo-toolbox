import sys
import string
import csv
import operator
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
from time import clock

sys.path.append('/bubo/home/h9/mircomic/glob/mypython/networkx-1.7/build/lib')
import networkx as nx
sys.path.append('/bubo/home/h9/mircomic/lib64/python')
import pygraphviz as pgv

import parse_dom_assign_hmmer3_repeat
import parse_hmmer3_domtblout
import parse_pdb_index
import parse_pdb_entry_type
import dict_pfam_pdb as pfam_pdb
#import dict_pfam_pdb_clan as pfam_pdb
import dict_pfam_numseq as pfam_numseq
import dict_pfam_numseq_clan as clan_numseq
#from domains import prot_dom_dict

#from pdb_seqres_2013 import scan_result
from pdb_seqres_2013_clan import scan_result


def contains(combi, target_combi):
    c_len = len(combi)
    tc_len = len(target_combi)
    len_diff = tc_len - c_len
    if len_diff < 0:
        return False
    for i in xrange(len_diff + 1):
        snippet = target_combi[i:(c_len + i)]
        if combi == snippet:
            return True
    return False


def create_dom_dict(domfile, pfam_clan, clan_name, pfameval, extraeval, pdb_res_dict):

    ## parse hmmscan result file into two dictionaries
    repeat_flag = True
    t0 = clock()
    #scan_result = parse_dom_assign_hmmer3_repeat.main(domfile, pfam_clan, clan_name, pfameval, extraeval)
    #scan_result = parse_hmmer3_domtblout.main(domfile, pfameval, extraeval, repeat_flag, pfam_clan, clan_name)
    #scan_result = parse_hmmer3_domtblout.main(domfile, pfameval, extraeval, repeat_flag)
    t1 = clock()
    print 'Parsing complete in %ds' % (t1 - t0)
    domfile.close()
    #prot_dom_dict = scan_result[0]
    prot_dom_dict = scan_result[2]
    dom_prot_dict = scan_result[1]
    dom_name_dict = scan_result[0]

    print len(prot_dom_dict)

    resolution_dict = parse_pdb_index.get_resolutions()
    entry_type_dict = parse_pdb_entry_type.get_entry_type()
    
    print '---'
    print len(resolution_dict)
    print len(entry_type_dict)
    print '---'
    
    prot_dict = {}
    prot_counts = {}
    single_dom_prot_dict = defaultdict(list)

    for prot_id, dom_list  in prot_dom_dict.iteritems():
        
        prot_counts[prot_id] = len(dom_list)
        
        pdb_acc = prot_id.split('_')[0]
        if not pdb_acc in resolution_dict:
            continue
        if resolution_dict[pdb_acc] <= 0 or resolution_dict[pdb_acc] > 2.5:
            continue
        if not pdb_acc in entry_type_dict:
            continue
        mol_type, meth_type = entry_type_dict[pdb_acc]
        if not mol_type == 'prot':
            continue
        
        if len(dom_list) == 1:
            single_dom_prot_dict[dom_list[0][2]].append((prot_id, resolution_dict[pdb_acc]))
            continue
        prot_dict[prot_id] = dom_list

    print len(prot_dict)
    #print prot_dict.keys()
    #print prot_dom_dict['1n11_A']
    
    prot_counts_sorted = sorted(prot_counts.iteritems(), key=operator.itemgetter(1), reverse=True)
    print prot_counts_sorted[:20]

    dom_counts = dict([(dom_id, len(seq_list)) for dom_id, seq_list in dom_prot_dict.iteritems()])
    #print max(dom_counts.iteritems(), key=operator.itemgetter(1))
    dom_counts_sorted = sorted(dom_counts.iteritems(), key=operator.itemgetter(1), reverse=True)
    print dom_counts_sorted[:20]

    dom_combi_dict = defaultdict(list)
    combi_prot_dict = defaultdict(list)
    combi_counts = defaultdict(int)
    combi_seqlen_dict = defaultdict(list)
    #for prot_id, dom_list in prot_dom_dict.iteritems():
    for prot_id, dom_list in prot_dict.iteritems():
        if len(dom_list) <= 1:
            continue
        combi = tuple([dom[2] for dom in dom_list])
        combi_counts[combi] += 1
        combi_pos = [dom[0] for dom in dom_list]
        combi_prot_dict[combi].append((prot_id, resolution_dict[prot_id.split('_')[0]], combi_pos))
        combi_seqlen = dom_list[-1][0][1] - dom_list[0][0][0]
        combi_seqlen_dict[combi].append(combi_seqlen)
        for dom in dom_list:
            dom_combi_dict[dom[2]].append(combi)

  
    combi2_dict = defaultdict(list)
    for combi in combi_counts.keys():
        for target_combi in combi_counts.keys():
            if len(combi) >= len(target_combi):
                continue
            if contains(combi, target_combi):
                combi2_dict[combi].append(target_combi)

    combi_counts_sorted = sorted(combi_counts.iteritems(), key=operator.itemgetter(1), reverse=True)
    #combifile = open('combi_counts_hq.csv', 'w')
    #combifile = open('combi_counts_hq_n4.csv', 'w')
    combifile = open('combi_counts_hq_l.csv', 'w')
    queryfile = open('query.txt', 'w')
    for (combi, count) in combi_counts_sorted:
        #combi_numseq = tuple([pfam_numseq.pfam_numseq_dict[dom_id] for dom_id in combi])
        combi_numseq = []
        for dom_id in combi:
            if dom_id in clan_numseq.pfam_numseq_dict:
                combi_numseq.append(clan_numseq.pfam_numseq_dict[dom_id])
            else:
                combi_numseq.append(pfam_numseq.pfam_numseq_dict[dom_id])
        has_enough_seq = min(combi_numseq) >= 500
        combi_names = tuple([dom_name_dict[dom_id] for dom_id in combi])
        combi_single_pdbs = tuple([dom_id in single_dom_prot_dict for dom_id in combi])
        #combi_pdbs = sorted(combi_prot_dict[combi], key=lambda x: x[1])
        combi_best_pdb = min(combi_prot_dict[combi], key=lambda x: x[1])
        combi2_names = []
        combi_seqlen = sum(combi_seqlen_dict[combi]) / float(len(combi_seqlen_dict[combi]))
        for combi2 in combi2_dict[combi]:
            combi2_names.append(tuple([dom_name_dict[dom_id] for dom_id in combi2]))
        #if all(combi_single_pdbs) and len(set(combi)) > 1 and len(combi2_dict[combi]) != 0 and has_enough_seq:
        #    combifile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (combi, combi_names, combi_single_pdbs, tuple(combi_numseq), count, len(combi2_dict[combi]), combi2_names, combi2_dict[combi], combi_pdbs))
        if all(combi_single_pdbs) and len(set(combi)) > 1 and has_enough_seq and len(combi) >= 3 and combi_seqlen > 50 and combi_seqlen < 500:
            combifile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (combi, combi_names, combi_single_pdbs, tuple(combi_numseq), count, combi_best_pdb, combi_seqlen))
            queryfile.write('%s\t%s\t%s\n' % (combi_best_pdb[0], combi_best_pdb[2], combi_names))
            #combifile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (combi, combi_names, combi_pdbs, count, len(combi2_dict[combi]), combi2_names, combi2_dict[combi]))
    combifile.close()
    queryfile.close()

    


    dom_combi_dict_unique = defaultdict(list)
    for dom_id, combi_list in dom_combi_dict.iteritems():
        dom_combi_dict_unique[dom_id] = list(set(combi_list))

    dom_combilen_dict = {}
    for dom_id, combi_list in dom_combi_dict.iteritems():
        avg_len = 0
        for combi in combi_list:
            avg_len += len(combi) / float(len(combi_list))
        dom_combilen_dict[dom_id] = avg_len

    dom_combilen_dict_unique = {}
    for dom_id, combi_list in dom_combi_dict_unique.iteritems():
        avg_len = 0
        for combi in combi_list:
            avg_len += len(combi) / float(len(combi_list))
        dom_combilen_dict_unique[dom_id] = avg_len


    dom_combi_counts = dict([(dom_id, len(combi_list)) for dom_id, combi_list in dom_combi_dict.iteritems()])
    dom_combi_counts_sorted = sorted(dom_combi_counts.iteritems(), key=operator.itemgetter(1), reverse=True)
    dom_combi_counts_unique = dict([(dom_id, len(combi_list)) for dom_id, combi_list in dom_combi_dict_unique.iteritems()])
    dom_combi_counts_unique_sorted = sorted(dom_combi_counts_unique.iteritems(), key=operator.itemgetter(1), reverse=True)
    print dom_combi_counts_sorted[:20]
    print dom_combi_counts_unique_sorted[:20]
    #print dom_combi_dict['PF00069.20']
    #print dom_combi_dict['PF07654.10'][0:10]

    unique_frac = defaultdict(list)
    for dom_id, count in dom_combi_counts.iteritems():
        frac = dom_combi_counts_unique[dom_id] / float(dom_combi_counts[dom_id])
        unique_frac[dom_id] = frac
    mean_unique_frac = sum(unique_frac.values()) / len(unique_frac)
    print mean_unique_frac

    #outfile = open('pdb_seqres_2013-02_2_multidomains.csv', 'w')
    outfile = open('tmp.csv', 'w')
    outfile.write('ID\tnum_unique_comb\tavg_unique_comb_len\tnum_comb\tavg_comb_len\tfraction\tsingle_struct\tName\n')
    for dom_id, count_unique in dom_combi_counts_unique_sorted:
        avg_combilen_unique = dom_combilen_dict_unique[dom_id]
        count = dom_combi_counts[dom_id]
        avg_combilen = dom_combilen_dict[dom_id]
        frac = float(count_unique) / float(count)
        num_single = len(single_dom_prot_dict[dom_id])
        name = dom_name_dict[dom_id]
        line = '%s\t%d\t%s\t%d\t%s\t%s\t%r\t%s\n' % (dom_id, count_unique, avg_combilen_unique, count, avg_combilen, frac, num_single, name)
        outfile.write(line)
    outfile.close()

    dom_idx_dict = dict([(entry[0], idx) for idx, entry in enumerate(dom_combi_counts_unique_sorted)])
    idx_dom_dict = dict([(idx, entry[0]) for idx, entry in enumerate(dom_combi_counts_unique_sorted)])

    partners = np.zeros((len(dom_combi_counts_unique_sorted), len(dom_combi_counts_unique_sorted)))
    #partners = np.zeros((100,100))
    for dom_id, combi_list in dom_combi_dict_unique.iteritems():
        for combi in combi_list:
            for dom_id2 in combi:
                idx1 = dom_idx_dict[dom_id]
                idx2 = dom_idx_dict[dom_id2]
                #if idx1 >= 100 or idx2 >= 100:
                #    continue
                if idx1 < idx2:
                    partners[idx1][idx2] += 1
                elif idx1 > idx2:
                    partners[idx2][idx1] += 1
                #if (idx1, idx2) not in G.edges() and (idx2, idx1) not in G.edges():
                #    G.add_edge(idx1, idx2)

    print partners
    np.savetxt('partners.csv', partners, delimiter=',')

    G = nx.Graph()
    for (i,j), value in np.ndenumerate(partners):
        if value != 0 and (i in range(10) and j in range(10)):# and i < 100 and j < 100:
            #print '%s - %s: %s' % (i, j, value)
            #G.add_edge(i,j,{'num_comb': value})
            G.add_edge(dom_name_dict[idx_dom_dict[i]],dom_name_dict[idx_dom_dict[j]])

    print len(G.nodes())
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #heatmap = ax.pcolor(partners, cmap=plt.cm.Blues)
    #pp = PdfPages('partners.pdf')
    #pp.savefig(fig)
    #pp.close()
    #plt.show()
    #pos = nx.shell_layout(G)
    pos = nx.spring_layout(G, scale=50)
    #pos = nx.random_layout(G)
    nx.draw(G, pos, font_size=4, node_size=0, node_color='#A0CBE2', edge_color='#BB0000', edge_width=1)#, width=1, edge_cmap=plt.cm.Blues)
    #A = nx.to_agraph(G)
    #A.layout()
    #A.draw('partners.pdf')
    #plt.show()
    plt.savefig('partners.pdf')#, dpi=1000) 

    """
    ## get domain statistics
    dom_dict = {}  # [len, avg E-val, avg Bit-score, max Bit-score, avg #rep, min #rep, max #rep, #seq, best PDB, longest PDB, name]
    
    for prot_id, dom_list  in prot_dom_dict.iteritems():
        dom_counts_dict = {}
        dom_pos_dict = {}
        dom_pdb = prot_id # since hmmscan was run on PDB sequences
                
        if len(dom_list) > 0:

            for dom in dom_list:
                dom_pos = dom[0]
                dom_id = dom[2]
                if dom_pos_dict.has_id(dom_id):
                    dom_pos_dict[dom_id].append(dom_pos)
                    dom_pos_dict[dom_id].sort()
                else:
                    dom_pos_dict[dom_id] = [dom_pos]
                if dom_counts_dict.has_id(dom_id):
                    dom_counts_dict[dom_id] += 1
                else:
                    dom_counts_dict[dom_id] = 1


            for dom in dom_list:
                if pdb_res_dict.has_id(dom_pdb[0:4]):
                    dom_pos = dom[0]
                    #dom_eval = dom[1] # not used here
                    dom_id = dom[2]
                    dom_len = dom[3]
                    dom_score = dom[4]
                    dom_count = dom_counts_dict[dom_id]
                    dom_res = pdb_res_dict[dom_pdb[0:4]]
                    
                    dom_clan = 'None'
                    if pfam_clan.has_id(string.split(dom_id,'.')[0]):
                        dom_clan = pfam_clan[string.split(dom_id,'.')[0]]
                    dom_name = domkey_name_dict[dom_id]

                    if dom_dict.has_id(dom_id):
                        dom_dict[dom_id][0] = max(dom_dict[dom_id][0], dom_len)
                        #dom_dict[dom_id][2] = (dom_dict[dom_id][2] + dom_score) / 2.0

                        # select best structure
                        old_res = pdb_res_dict[dom_dict[dom_id][6][0:4]]
                        old_score = dom_dict[dom_id][2]
                        if dom_score > old_score:
                            dom_dict[dom_id][2] = dom_score 
                        if dom_res != -1 and old_res != -1:
                            #if (dom_score / dom_res) < (old_score / old_res):
                            if dom_res < 2.8 and dom_score > old_score:
                                dom_dict[dom_id][1] = dom_score 
                                dom_dict[dom_id][6] = dom_pdb
                                dom_dict[dom_id][7] = dom_res
                                dom_dict[dom_id][8] = '%d-%d' % (dom_pos[0], dom_pos[1])
                        elif dom_res != -1:
                            dom_dict[dom_id][1] = dom_score 
                            dom_dict[dom_id][6] = dom_pdb
                            dom_dict[dom_id][7] = dom_res
                            dom_dict[dom_id][8] = '%d-%d' % (dom_pos[0], dom_pos[1])
                        elif dom_res == -1 and old_res == -1:
                            if dom_score > old_score:
                                dom_dict[dom_id][1] = dom_score 
                                dom_dict[dom_id][6] = dom_pdb
                                dom_dict[dom_id][7] = dom_res
                                dom_dict[dom_id][8] = '%d-%d' % (dom_pos[0], dom_pos[1])

                        # select structure with most repeats
                        old_max_count = dom_dict[dom_id][5]
                        #old_min_pos = int(dom_dict[dom_id][11].split('-')[0])
                        #old_max_pos = int(dom_dict[dom_id][11].split('-')[1])
                        dom_dict[dom_id][3] = (dom_dict[dom_id][3] + dom_count) / 2
                        dom_dict[dom_id][4] = min(dom_dict[dom_id][4], dom_count)
                        dom_dict[dom_id][5] = max(old_max_count, dom_count)
                        if dom_count > old_max_count:
                            dom_dict[dom_id][9] = dom_pdb
                            dom_dict[dom_id][10] = dom_res
                            dom_dict[dom_id][11] = dom_pos_dict[dom_id]
                            #dom_dict[dom_id][11] = '%d-%d' % (min(old_min_pos, dom_pos[0]),
                            #                                   max(old_max_pos, dom_pos[1]))
                        elif dom_count == old_max_count:
                            if dom_res < old_res:
                                dom_dict[dom_id][9] = dom_pdb
                                dom_dict[dom_id][10] = dom_res
                                dom_dict[dom_id][11] = dom_pos_dict[dom_id]
                                #dom_dict[dom_id][11] = '%d-%d' % (min(old_min_pos, dom_pos[0]),
                                #                                   max(old_max_pos, dom_pos[1]))


                    else:
                        dom_dict[dom_id] = [dom_len,
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
                                            dom_pos_dict[dom_id],
                                            dom_clan,
                                            dom_name]
        """        


def write_dom_dict(dom_dict, outfile):
    #sorted(pfam_stat_dict.items(), key = lambda x : x[1][1], reverse = True)[0:10]

    writer = csv.writer(outfile, delimiter='\t')
    colnames = ['PfamID', 'length', 'score_in_best_PDB', 'max_score', 'avg_num_repeats', 'min_num_repeats', 'max_num_repeats', 'best_PDB', 'res_in_best_PDB', 'pos_in_best_PDB', 'longest_PDB', 'res_in_longest_PDB', 'pos_in_longest_PDB', 'ClanID', 'Name']
    writer.writerow(colnames)

    for dom_id, dom_statlist in dom_dict.iteritems():
        line = [dom_id] + dom_statlist
        #print [dom_id] + dom_statlist
        #print line
        writer.writerow(line)


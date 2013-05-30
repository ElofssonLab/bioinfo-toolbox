import sys
import string
import csv

#import parse_dom_assign_hmmer3_repeat
import parse_hmmer3_domtblout
#import dict_pfam_pdb as pfam_pdb
import dict_pfam_pdb_clan as pfam_pdb
#import dict_pfam_numseq as pfam_numseq
import dict_pfam_numseq_clan as pfam_numseq

def create_dom_dict(domfile, pfam_clan, clan_name, pfameval, extraeval):

    ## parse hmmscan result file into two dictionaries
    #scan_result = parse_dom_assign_hmmer3_repeat.main(domfile, pfam_clan, clan_name, pfameval, extraeval)
    repeat_flag = True
    scan_result = parse_hmmer3_domtblout.main(domfile, pfameval, extraeval, repeat_flag, pfam_clan, clan_name)
    domfile.close()

    prot_dom_dict = scan_result[2]
    domkey_name_dict = scan_result[0]
    
    ## get domain statistics
    #dom_dict = {}  # [len, avg E-val, avg #rep, min #rep, max #rep, #seq, PDB, clan]
    dom_dict = {}  # [len, avg E-val, avg Bit-score, avg #rep, min #rep, max #rep, #seq, PDB, name]
    
    for prot_key, dom_list  in prot_dom_dict.iteritems():
        dom_counts_dict = {}
                
        if len(dom_list) > 0:

            for dom in dom_list:
                dom_pos = dom[0]
                dom_len = abs(dom_pos[1] - dom_pos[0])
                #dom_len = dom[3]
                dom_eval = dom[1]
                dom_score = dom[3]
                dom_key = dom[2]
                if pfam_numseq.pfam_numseq_dict.has_key(dom_key):
                    dom_numseq = pfam_numseq.pfam_numseq_dict[dom_key]
                else:
                    continue
                
                #dom_pdb = int(pfam_pdb.pfam_pdb_dict.has_key(dom_key))
                dom_pdb = 'None'
                if pfam_pdb.pfam_pdb_dict.has_key(dom_key):
                    dom_pdb = pfam_pdb.pfam_pdb_dict[dom_key]
                
                dom_clan = 'None'
                #if pfam_clan.has_key(string.split(dom_key,'.')[0]):
                #    dom_clan = pfam_clan[string.split(dom_key,'.')[0]]
                #dom_name = domkey_name_dict[dom_key]
                if dom_key.startswith('CL'):
                    dom_clan = dom_key
                dom_name = domkey_name_dict[dom_key]

                if dom_dict.has_key(dom_key):
                    dom_dict[dom_key][0] = max(dom_dict[dom_key][0], dom_len)
                    dom_dict[dom_key][1] = (dom_dict[dom_key][1] + dom_eval) / 2.0
                    dom_dict[dom_key][2] = (dom_dict[dom_key][2] + dom_score) / 2.0
                else:
                    dom_dict[dom_key] = [dom_len,
                                        dom_eval,
                                        dom_score,
                                        0.0,
                                        sys.maxint,
                                        0,
                                        dom_numseq,
                                        dom_pdb,
                                        dom_clan,
                                        dom_name
                                        ]
        
                if dom_counts_dict.has_key(dom_key):
                    dom_counts_dict[dom_key] += 1
                else:
                    dom_counts_dict[dom_key] = 1

            for dom_counts_key, dom_count in dom_counts_dict.iteritems():
                if dom_dict[dom_counts_key][3] != 0.0:
                    dom_dict[dom_counts_key][3] = (dom_dict[dom_counts_key][3] + dom_count) / 2.0
                else:
                    dom_dict[dom_counts_key][3] = dom_count
                
                dom_dict[dom_counts_key][4] = min(dom_dict[dom_counts_key][4], dom_count) 
                dom_dict[dom_counts_key][5] = max(dom_dict[dom_counts_key][5], dom_count) 

    return dom_dict


def write_dom_dict(dom_dict, outfile):
    #sorted(pfam_stat_dict.items(), key = lambda x : x[1][1], reverse = True)[0:10]

    writer = csv.writer(outfile, delimiter='\t')
    colnames = ['PfamID', 'length', 'avg_E-value', 'avg_score', 'avg_num_repeats', 'min_num_repeats', 'max_num_repeats', 'num_sequences', 'PDB', 'ClanID', 'Name']
    writer.writerow(colnames)

    for dom_key, dom_statlist in dom_dict.iteritems():
        line = [dom_key] + dom_statlist
        #print [dom_key] + dom_statlist
        #print line
        writer.writerow(line)


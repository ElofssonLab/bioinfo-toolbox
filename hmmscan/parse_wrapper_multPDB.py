import sys

import analyse_hmmscan_result as analyse
import dict_pfam_clan as pfam_clan
import dict_clan_name as clan_name

#import dict_pfam_stat as pfam_stat
import parse_pdb_index


domfile = open(sys.argv[1], 'r')
#pfam_clan = {'PF00880': 'PF00880'}
#clan_name = {'PF00880':'Nebulin'}
#pfameval = 10   ##had set at 0.01
#extraeval = 10              ##more permissive for predictions of IG-like fildom
pfameval = 0.1   ##had set at 0.01
extraeval = 1.0              ##more permissive for predictions of IG-like fildom
#lista = parse_dom_assign_hmmer3_repeat.main(domfile, pfam_clan, clan_name, pfameval, extraeval)

rep_cutoff = 3

pdb_resolu_dict = parse_pdb_index.get_resolutions()
dom_dict = analyse.create_dom_dict(domfile, pfam_clan.pfam_clan_dict, clan_name.clan_name_dict, pfameval, extraeval, pdb_resolu_dict)
#dom_dict = pfam_stat.pfam_stat_dict
domfile.close()

#outfile = open(sys.argv[2], 'wb')
#analyse.write_dom_dict(dom_dict, outfile)
#outfile.close()

#print dom_dict
#print len(dom_dict)
#print len(lista[0])
#print lista[0]

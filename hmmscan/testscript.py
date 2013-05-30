import sys

#import parse_pfam_full as parse
import dict_pfam_clan as pfam_clan

#import parse_pfam_clan as parse
import parse_pfam_pdb as parse

domfile = open(sys.argv[1], 'r')
#lista = parse.main_clan(domfile, pfam_clan.pfam_clan_dict)
lista = parse.main(domfile)
#lista,listb = parse.main(domfile)
domfile.close()

print lista
#print len(lista)

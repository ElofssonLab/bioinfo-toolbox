import sys
import string, copy

## parses Pfam-C
## output: dictionary {pfam_acc: clan_acc}
def main(afile):

    pfam_clan_dict = {}
    clan_name_dict = {}

    aline = afile.readline()

    clan_acc = ''
    clan_name = ''

    while not aline == '':
        alist = string.split(aline)

        if alist[0] == '#=GF':
            if alist[1] == 'ID' and clan_acc != '':
                clan_name = string.strip(alist[2])
                clan_name_dict[clan_acc] = clan_name
            elif alist[1] == 'ID' and clan_acc == '':
                clan_name = string.strip(alist[2])
            if alist[1] == 'AC' and clan_name != '':
                clan_acc = string.strip(alist[2])
                clan_name_dict[clan_acc] = clan_name
            elif alist[1] == 'AC' and clan_name == '':
                clan_acc = string.strip(alist[2])
            if alist[1] == 'MB' and clan_acc != '':
                pfam_acc = string.strip(alist[2],';')
                pfam_clan_dict[pfam_acc] = clan_acc
        if alist[0] == '//':
            clan_acc = ''
            clan_name = ''

        aline = afile.readline()

    return pfam_clan_dict, clan_name_dict


if __name__ == '__main__':

    pfam_clan_dict, clan_name_dict = main(open(sys.argv[1], 'r'))
    print clan_name_dict

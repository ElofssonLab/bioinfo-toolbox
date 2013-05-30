import string, copy

## parses 'pdb_pfam_mapping.txt'
## output: dictionary {pfam_acc: pdb_acc}
def main(afile):

    pfam_pdb_dict = {}

    aline = afile.readline()  ##skip file header
    aline = afile.readline()

    while not aline == '':
        alist = string.split(aline)

        pfam_acc = alist[4]
        pdb_acc = alist[0]
        pfam_pdb_dict[pfam_acc] = pdb_acc

        aline = afile.readline()

    return pfam_pdb_dict

## parses 'pdb_pfam_mapping.txt' on clan level
## output: dictionary {pfam_acc/clan_acc: pdb_acc}
def main_clan(afile, pfam_clan):

    pfam_pdb_dict = {}

    aline = afile.readline()  ##skip file header
    aline = afile.readline()

    while not aline == '':
        alist = string.split(aline)
        
        pfam_acc = alist[4]

        if pfam_clan.has_key(string.split(pfam_acc,'.')[0]):
            pfam_acc = pfam_clan[string.split(pfam_acc,'.')[0]]
        
        pdb_acc = alist[0]
        pfam_pdb_dict[pfam_acc] = pdb_acc

        aline = afile.readline()

    return pfam_pdb_dict


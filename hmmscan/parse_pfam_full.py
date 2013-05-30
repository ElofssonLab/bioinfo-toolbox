import string, copy

## parses 'Pfam-A.full'
## output: dictionary {pfam_acc: number_of_sequences}
def main(afile):

    pfam_numseq_dict = {} 

    aline = afile.readline()

    pfam_acc = ''

    while not aline == '':
        alist = string.split(aline)

        if alist[0] == '#=GF':
            if alist[1] == 'AC':
                pfam_acc = string.strip(alist[2])
                #if pfam_clan.has_key(string.split(pfam_acc,'.')[0]):
                #    pfam_acc = pfam_clan[string.split(pfam_acc,'.')[0]]
            if alist[1] == 'SQ' and pfam_acc != '':
                pfam_numseq_dict[pfam_acc] = string.atoi(string.strip(alist[2]))

        aline = afile.readline()

    return pfam_numseq_dict


## parses 'Pfam-A.full' on clan level
## output: dictionary {pfam_acc|clan_acc: number_of_sequences}
def main_clan(afile, pfam_clan):

    pfam_numseq_dict = {} 

    aline = afile.readline()

    pfam_acc = ''

    while not aline == '':
        alist = string.split(aline)

        if alist[0] == '#=GF':
            if alist[1] == 'AC':
                pfam_acc = string.strip(alist[2])
                if pfam_clan.has_key(string.split(pfam_acc,'.')[0]):
                    pfam_acc = pfam_clan[string.split(pfam_acc,'.')[0]]
            if alist[1] == 'SQ' and pfam_acc != '':
                if pfam_numseq_dict.has_key(pfam_acc):
                    pfam_numseq_dict[pfam_acc] += string.atoi(string.strip(alist[2]))
                else:
                    pfam_numseq_dict[pfam_acc] = string.atoi(string.strip(alist[2]))


        aline = afile.readline()

    return pfam_numseq_dict


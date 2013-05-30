import sys
import urllib

sys.path.append('/bubo/sw/apps/bioinfo/biopython/1.56/lib/python')

from Bio import SeqIO
from Bio import ExPASy
from Bio import SwissProt

import parse_fasta
import check_sequence_overlap as check


def overlap_filter(afile, query_id):

    ali_dict = parse_fasta.read_fasta(afile, query_id)

    query_seq = ali_dict[query_id][0]

    counts = {0:0, 1:0, 2:0, 3:0, 4:0}

    for header, seqlist in ali_dict.iteritems():
        nseq = len(seqlist)
        counts[3] += nseq

        if nseq > 1 and header != query_id:
            #print nseq
            target_acc = header.split('|')[1]

            try:
                handle = urllib.urlopen('http://www.uniprot.org/uniprot/%s.txt' % target_acc)
                target_record = SwissProt.read(handle)
                target_seq = target_record.sequence
                #print target_acc
                for i in range(nseq):
                    for j in range(i + 1, nseq):
                        #flag = check.check_jackhmmer(seqlist[i], seqlist[j],
                        #                             target_acc, target_seq, query_seq)
                        flag = check.check_fast(seqlist[i], seqlist[j],
                                                target_acc, target_seq, query_seq)
                        counts[flag] += 1
                        if flag == 1:
                            del ali_dict[header][1]
                        if flag == 2:
                            del ali_dict[header][0]
            except ValueError:
                counts[4] += 1

    return counts, ali_dict


def overlap_filter_jackhmmer(afile, query_id):

    ali_dict = parse_fasta.read_fasta(afile, query_id)
    ali_dict_jhmmer = {}
    output_dict = {}

    query_seq = ali_dict[query_id][0]

    counts = {0:0, 1:0, 2:0, 3:0, 4:0}

    for header, seqlist in ali_dict.iteritems():
        if header == query_id:
            pos = (1, len(seqlist[0]))
            ali_dict_jhmmer[header] = {pos : seqlist[0]}
        else:
            #acc = header.split('|')[1]
            pos_string = header.split('|')[2].split()[0].split('/')[-1]
            acc = header.replace(pos_string, '<POSITION>')
            pos = (int(pos_string.split('-')[0]), int(pos_string.split('-')[1]))
            if ali_dict_jhmmer.has_key(acc):
                ali_dict_jhmmer[acc][pos] = seqlist[0]
            else:
                ali_dict_jhmmer[acc] = {pos : seqlist[0]}

    for acc, pos_dict in ali_dict_jhmmer.iteritems():

        nseq = len(pos_dict)
        counts[3] += nseq
        #print nseq

        if nseq > 1 and acc != query_id:
            nonoverlap_list, overlap_count = check.resolve_overlaps(pos_dict, query_seq)
            
            new_pos_dict = {}
            for pos, seq in pos_dict.iteritems():
                if pos in nonoverlap_list:
                    new_pos_dict[pos] = seq
            output_dict[acc] = new_pos_dict

            counts[0] += len(nonoverlap_list)
            counts[1] += overlap_count

        elif acc == query_id:
            print acc
            output_dict[acc] = ali_dict_jhmmer[acc]
            print output_dict[acc]

    return counts, output_dict


def coverage_filter(afile, query_id, coverage):

    ali_dict = parse_fasta.read_fasta(afile, query_id)
    ref_len = len(ali_dict[query_id][0].translate(None, '.-').upper())

    output_dict = {}
    count = 0

    for header, seqlist in ali_dict.iteritems():
        new_seqlist = []
        for seq in seqlist:
            seq_len = len(seq.translate(None, '.-').upper())
            len_ratio = seq_len / float(ref_len)
            if len_ratio >= coverage:
                new_seqlist.append(seq)
            else:
                count += 1
        output_dict[header] = new_seqlist

    return count, output_dict



if __name__ == "__main__":
    args = sys.argv

    jflag = False
    cflag = False
    coverage = 1.0

    if '-j' in args:
        jflag = True
        args.remove('-j')

    if '-c' in args:
        cflag = True
        coverage_i = args.index('-c') + 1
        coverage = float(args[coverage_i])
        del args[coverage_i]
        args.remove('-c')

    afile = open(args[1], 'r')
    query_id = args[2]
    #target_db = SeqIO.index('~/databases/uniprot_sp_tr.fa', 'fasta')
    #counts, ali_dict = main(afile, query_id)

    if jflag and not cflag:
        counts, ali_dict = overlap_filter_jackhmmer(afile, query_id)
        afile.close()

        print counts

        outfile = open(args[3], 'w')
        
        # write query entry
        outfile.write('>' + query_id + '\n' + ali_dict[query_id].itervalues().next() + '\n')

        for header, pos_dict in ali_dict.iteritems():
            for pos, seq in pos_dict.iteritems():
                if header != query_id:
                    header_line = '>%s' % header
                    pos_string = '%d-%d' % pos
                    outfile.write(header_line.replace('<POSITION>', pos_string) + '\n' + seq + '\n')

    elif not jflag and not cflag:
        counts, ali_dict = overlap_filter(afile, query_id)
        afile.close()

        print counts

        outfile = open(args[3], 'w')

        # write query entry
        outfile.write('>' + query_id + '\n' + ali_dict[query_id][0] + '\n')

        for header, seqlist in ali_dict.iteritems():
            for seq in seqlist:
                if header != query_id:
                    header_line = '>%s' % header
                    outfile.write(header_line + '\n' + seq + '\n')

    elif cflag:
        count, ali_dict = coverage_filter(afile, query_id, coverage)
        afile.close()

        print count

        outfile = open(args[3], 'w')

        # write query entry
        outfile.write('>' + query_id + '\n' + ali_dict[query_id][0] + '\n')

        for header, seqlist in ali_dict.iteritems():
            for seq in seqlist:
                if header != query_id:
                    header_line = '>%s' % header
                    outfile.write(header_line + '\n' + seq + '\n')




    





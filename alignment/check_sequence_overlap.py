import subprocess
import sys
import urllib

sys.path.append('/bubo/sw/apps/bioinfo/biopython/1.56/lib/python')
from Bio import ExPASy
from Bio import SwissProt
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

def get_best_pos(overlap_list, pos_dict, query_seq):

    best_score = float('-inf')
    best_pos = overlap_list[0]

    for pos in overlap_list:
        query_seq = query_seq.translate(None, '.-').upper()
        dist_mat = MatrixInfo.blosum62

        seq = pos_dict[pos].translate(None, '.-').upper()
        score = pairwise2.align.globalds(query_seq, seq, dist_mat,
                                         -10, -0.5, score_only=True)
        if score > best_score:
            best_score = score
            best_pos = pos

    return best_pos



def resolve_overlaps(pos_dict, query_seq):

    nonoverlap_list = []
    overlap_count = 0

    pos_list = pos_dict.keys()
    pos_list.sort()
    #print pos_list
    while pos_list:
        curr_pos = pos_list[0]
        curr_overlaps = [pos for pos in pos_list
                         if pos[0] <= curr_pos[1]]
        #print curr_overlaps
        if len(curr_overlaps) > 1:
            overlap_count += len(curr_overlaps)
            best_pos = get_best_pos(curr_overlaps, pos_dict, query_seq)
            pos_list = [pos for pos in pos_list
                        if pos not in curr_overlaps or pos == best_pos]
        else:
            nonoverlap_list.append(curr_pos)
            pos_list.remove(curr_pos)

    #print nonoverlap_list
    return nonoverlap_list, overlap_count


# Checks whether two subsequences of the same UniProtKB entry are overlapping.
# If there is an overlap, the sequence aligning best to query_seq is indicated.
# output: 0 - no overlap; 1 - prefer seq1; 2 - prefer seq2
def check(seq1, seq2, target_acc, target_seq, query_seq):

    seq1 = seq1.translate(None, '.-').upper()
    seq2 = seq2.translate(None, '.-').upper()

    if seq1 == '' and not seq2 == '':
        result = 2
        return result
    elif not seq1 == '' and seq2 == '':
        result = 1
        return result
    elif seq1 == '' and seq2 == '':
        result = 4
        return result


    #print target_acc
    #handle = ExPASy.get_sprot_raw(target_acc)
    #print 'http://www.uniprot.org/uniprot/%s' % target_acc
    #handle = urllib.urlopen('http://www.uniprot.org/uniprot/%s.txt' % target_acc)
    #handle = urllib.urlopen('http://www.uniprot.org/uniprot/F6ZI51.txt')
    #handle = urllib.urlopen('http://www.uniprot.org/uniprot/D3ZEN4.txt')
    try:
        #target_record = SwissProt.read(handle)
        #target_seq = target_record.sequence

        #target_seq = db_dict[target_acc].seq
        
        #print 'target_seq=' + target_seq
        #print 'seq1=' + seq1
        #print 'seq2=' + seq2

        align1 = pairwise2.align.globalxs(target_seq, seq1, -10, -0.5)
        align2 = pairwise2.align.globalxs(target_seq, seq2, -10, -0.5)
        
        #print align1
        #print align2

        seq1_align = align1[0][1]
        seq2_align = align2[0][1]
        
        is_overlap = False

        for i in range(len(seq1_align)):
            if (seq1_align[i] != '-' and seq2_align[i] != '-'):
                is_overlap = True
                break

        result = 0

        if is_overlap:
            query_seq = query_seq.translate(None, '.-').upper()
            dist_mat = MatrixInfo.blosum62

            score1 = pairwise2.align.globalds(query_seq, seq1, dist_mat,
                                                    -10, -0.5, score_only=True)
            score2 = pairwise2.align.globalds(query_seq, seq2, dist_mat,
                                                    -10, -0.5, score_only=True)
            if score1 > score2:
                result = 1
            else:
                result = 2

    except ValueError:
        result = 4

    return result

    #print seq1_align
    #print seq2_align


def check_fast(seq1, seq2, target_acc, target_seq, query_seq):
    seq1 = seq1.translate(None, '.-').upper()
    seq2 = seq2.translate(None, '.-').upper()

    if seq1 == '' and not seq2 == '':
        result = 2
        return result
    elif not seq1 == '' and seq2 == '':
        result = 1
        return result
    elif seq1 == '' and seq2 == '':
        result = 4
        return result

    pos1x = target_seq.find(seq1)
    pos1y = pos1x + len(seq1)
    pos1 = [pos1x, pos1y]

    pos2x = target_seq.find(seq2)
    pos2y = pos2x + len(seq2)
    pos2 = [pos2x, pos2y]

    overlap = max(0, min(pos1[1], pos2[1]) - max(pos1[0], pos2[0]))

    if overlap != 0:
        query_seq = query_seq.translate(None, '.-').upper()
        dist_mat = MatrixInfo.blosum62

        score1 = pairwise2.align.globalds(query_seq, seq1, dist_mat,
                                                -10, -0.5, score_only=True)
        score2 = pairwise2.align.globalds(query_seq, seq2, dist_mat,
                                                -10, -0.5, score_only=True)
        if score1 > score2:
            result = 1
        else:
            result = 2
    else:
        result = 0

    return result
 



if __name__ == "__main__":
    #seq1 = '---VDDE---FTLIYASKGGYLELVNLLI----KNGADIHVNDD---APLKWASKNGHLEVVKYLVENGADIHAYNE----LVVYASEGGHLQIVKYLVKKGA----DIHAEDDE---ALKWASRSGHLEVVKYLVEKGANFR----AENDYALRWACEKGHLEIVKYLVEKGADIHAED---EYALRWASRSGHLEVVKYLVENGADIHACNDYG---LRKASRNRHLNVV---------'
    seq1 = '---VDDE---FTLIYASKGGYLELVNLLI----KNGADIHV'
    seq2 = '----NDD---APLKWASKNGHLEVVKYLVE----NGADIHAYNE----LVVYASEGGHLQIVKYLVKKGADIHAEDDE---ALKWASRSGHLEVVKYLVEKGA----NFR---AENDYALRWACEKGHLEIVKYLVEKGADIHAED----EYALRWASRSGHLEVVKYLVENGADIHACNDYG---LRKASRNRHLNVVKYLMENGANIHAKDDY---ALRLAS-----------------'
    acc = 'E3VXP1'
    target_seq = 'MISIGVYIHTRSEYPLRWASTEGHVEVVIYLVENGADLNVLKLFYVDKYLQVVKYLIENGSNIHVDDEFTLIYASKGGYLELVNLLIKNGADIHVNDDAPLKWASKNGHLEVVKYLVENGADIHAYNELVVYASEGGHLQIVKYLVKKGADIHAEDDEALKWASRSGHLEVVKYLVEKGANFRAENDYALRWACEKGHLEIVKYLVEKGADIHAEDEYALRWASRSGHLEVVKYLVENGADIHACNDYGLRKASRNRHLNVVKYLMENGANIHAKDDYALRLASVNGHFKLVKFLVENGANIHAKNNEALIRPLGEGHIKIVTYLE'
    query_seq = 'MATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPAS'
    check_result = check_fast(seq1, seq2, acc, target_seq, query_seq)

    print check_result


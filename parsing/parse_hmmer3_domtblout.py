import sys
import copy
from collections import defaultdict
from time import clock


def read(infile, dom_clan_dict={}, clan_name_dict={}):

    dom_name_dict = defaultdict(list)
    query_dom_dict = defaultdict(list)

    # remove first three commented lines
    #infile.readline()
    #infile.readline()
    #infile.readline()

    for line in infile:
        if line.startswith('#'):
            continue
        line_list = line.split()
        dom_name = line_list[0]
        dom_id = line_list[1]
        # when parsing non-pfam searches (e.g. pdb)
        if dom_id == '-':
            dom_id = dom_name
        dom_id_short = dom_id.split('.')[0]
        query_id = line_list[3]
        eval = float(line_list[12])
        bitscore = float(line_list[13])
        start = int(line_list[19])
        end = int(line_list[20])
        start_dom = int(line_list[15])
        end_dom = int(line_list[16])
        start_seq = int(line_list[17])
        end_seq = int(line_list[18])

        # if clan dictionaries are given, translate id and name
        try:
            clan_id = dom_clan_dict[dom_id_short]
            dom_name = clan_name_dict[clan_id]
            dom_id = clan_id
        except KeyError:
            pass

        query_dom_entry = ((start, end), eval, dom_id, bitscore,
                (start_dom, end_dom), (start_seq, end_seq))
        dom_name_dict[dom_id] = dom_name
        query_dom_dict[query_id].append(query_dom_entry)

    return [dict(dom_name_dict), dict(query_dom_dict)]


    
def remove_overlaps(query_dom_dict, repeat_flag):

    # datatype:
    # -> query_dom_dict = {query_id: entries}
    # -> entries = [((start, end), eval, dom_id, bitscore)]

    nonoverlap_dict = {}

    for query_id, entries in query_dom_dict.iteritems():
        # sort by start position of domain in the sequence
        entries.sort(key=lambda entry: entry[0][0])
        del_entries = []
        nonoverlap_entries = []
        while entries:
            curr_entry = entries[0]
            if repeat_flag:
                # get all entries that overlap with current entry
                # and have the same id (ONLY CHECK FOR REPEAT DOMAINS)
                curr_overlaps = [entry for entry in entries
                                 if entry[0][0] <= curr_entry[0][1]
                                 and entry[2] == curr_entry[2]]
            else:
                # get all entries that overlap with current entry
                curr_overlaps = [entry for entry in entries
                                 if entry[0][0] <= curr_entry[0][1]]
            # if there are overlaps:
            if len(curr_overlaps) > 1:
                # get entry with lowest evalue from all overlapping entries
                best_entry = min(curr_overlaps, key=lambda entry: entry[1])
                if best_entry != curr_entry:
                    # mark current entry for deletion, 
                    # keep the rest for further checking
                    del_entries = [curr_entry]
                else:
                    # keep current entry and mark all other overlapping entries
                    curr_overlaps.remove(curr_entry)
                    del_entries = curr_overlaps
                # delete entries according to above
                entries = [entry for entry in entries
                           if entry not in del_entries]
            else:
                # no overlapping entries found: keep current entry
                nonoverlap_entries.append(curr_entry)
                entries.remove(curr_entry)
        nonoverlap_dict[query_id] = nonoverlap_entries

    return nonoverlap_dict
       


def has_valid_neighbor(perm_entry, perm_idx, entries, strict_eval):
    
    pre_neighbors = entries[:perm_idx][::-1]
    post_neighbors = entries[(perm_idx + 1):]

    for pre in pre_neighbors:
        if pre[2] != perm_entry[2]:
            continue
        dist = abs(perm_entry[0][0] - pre[0][1])
        if dist > 50:
            break
        if pre[1] < strict_eval:
            return True

    for post in post_neighbors:
        if post[2] != perm_entry[2]:
            continue
        dist = abs(post[0][0] - perm_entry[0][1])
        if dist > 50:
            break
        if post[1] < strict_eval:
            return True

    return False



def evalue_filter(query_dom_dict, strict_eval, perm_eval):
 
    # datatype:
    # -> query_dom_dict = {query_id: entries}
    # -> entries = [((start, end), eval, dom_id, bitscore)]

    filtered_dict = {}
    rep_count = 0
    perm_rep_count = 0

    for query_id, entries in query_dom_dict.iteritems():
        # sort by start position of domain in the sequence
        entries.sort(key=lambda entry: entry[0][0])

        # remove domains with higher evalue than permissive
        entries = [entry for entry in entries
                   if entry[1] < perm_eval]
        # if there are no domains left -> ignore whole query sequence
        if not entries:
            continue

        # if there are no repeated domains (all domain ids are unique)
        # -> apply strict evalue filter to all domains
        dom_ids = map(lambda x: x[2], entries)
        unique_ids = list(set(dom_ids))
        strict_entries = [entry for entry in entries
                          if entry[1] < strict_eval]
        if not strict_entries:
            continue
        if len(dom_ids) == len(unique_ids):
            filtered_dict[query_id] = strict_entries
            continue
        
        # at this point 'entries' contains at least one duplicated domain
        rep_count += 1

        if len(strict_entries) == len(entries):
            filtered_dict[query_id] = strict_entries
            continue

        # at this point 'entries' also contains at least one domain
        # with an evalue in the permissive range
        perm_rep_count += 1
        final_entries = strict_entries
        
        perm_entries = dict([(idx, entry) for idx, entry in enumerate(entries)
                        if entry[1] >= strict_eval])
        
        # only if a permissive entry has a neighbor which has:
        # 1) a strict evalue AND
        # 2) is closer than 50 residues in the sequence
        # -> it is not deleted
        for perm_idx, perm_entry in perm_entries.iteritems():
            if has_valid_neighbor(perm_entry, perm_idx, entries, strict_eval):
                final_entries.append(perm_entry)

        final_entries.sort(key=lambda entry: entry[0][0])
        filtered_dict[query_id] = final_entries

    print 'Detected: - %d sequences with at least one duplicated domain.' % rep_count
    print '          - %d of which have at least one domain between your evalue thresholds.' % perm_rep_count

    return filtered_dict




def main(infile, strict_eval, perm_eval, repeat_flag, dom_clan_dict={}, clan_name_dict={}):

    print 'Parsing started...'
    t0 = clock()
    repeat_flag = False
    dom_query_dict = defaultdict(list)
    query_dom_dict = defaultdict(list)
    domtblout = read(infile, dom_clan_dict, clan_name_dict)
    t1 = clock()
    print 'Parsing ended in %ds.\n' % (t1 - t0)
    
    dom_name_dict = domtblout[0]
    raw_query_dom_dict = domtblout[1]

    print 'Removing overlaps...'
    t0 = clock()
    nonoverlap_query_dom_dict = remove_overlaps(raw_query_dom_dict, repeat_flag)
    t1 = clock()
    print 'Overlaps removed in %ds.\n' % (t1 - t0)

    print 'Running evalue filter...'
    t0 = clock()
    query_dom_dict = evalue_filter(nonoverlap_query_dom_dict, strict_eval, perm_eval)
    t1 = clock()
    print 'Filtered in %ds.\n' % (t1 - t0)

    for query_id, entries in query_dom_dict.iteritems():
        for entry in entries:
            dom_id = entry[2]
            dom_query_entry = (entry[0], entry[1], query_id, entry[3])
            dom_query_dict[dom_id].append(dom_query_entry)

    print 'There are %d different domains in %d sequences.' % (len(dom_query_dict), len(query_dom_dict))

    return [dict(dom_name_dict), dict(dom_query_dict), dict(query_dom_dict)]




if __name__ == "__main__":
    infile = open(sys.argv[1], 'r')
    domtblout = main(infile, 0.01, 0.1, False)
    infile.close()
    #print domtblout[0]
    #print domtblout[1]
    #print domtblout[2]
    print 'DONE!'
    
    if len(sys.argv) == 3:
        with open(sys.argv[2]) as ids:
            for id in ids:
                id = id.strip()
                if id in domtblout[2].keys():
                    minpos = min(domtblout[2][id], key=lambda x: x[1])
                    print id, minpos[0], minpos[1]
                else:
                    print id



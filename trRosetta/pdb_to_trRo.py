#!/usr/bin/env python3
import sys
import os
import numpy as np
import scipy.stats as st
from collections import defaultdict


def parse_atm_record(line):
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    
    return record


def get_first_chain(pdbfile):
    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        break

    return atm_record['chain']


def get_res_dict(pdbfile, chain):
    cb_lst = []
    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue

        atm_record = parse_atm_record(line)

        if atm_record['chain'] != ' ' and atm_record['chain'] != chain  and chain != '*':
            continue

        res_i = atm_record['res_no']
        
        if res_dict.keys():
            min_res_i = min(res_dict.keys())
        else:
            min_res_i = res_i

        # I do not really understand whatis tested here, but seems to never be true/AE
        if res_i > 1000 and len(res_dict) < 1000 and min_res_i + len(res_dict) < 1000:
            continue

        # Why ?
        if atm_record['insert'] == 'X':
            res_i = res_i * 0.001

        atm = [float('inf'), float('inf'), float('inf')]

        if atm_record['atm_name'] == 'CA':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm))   
        elif atm_record['atm_name'] == 'CB':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm)) 
    
    return res_dict


def get_cb_coordinates(pdbfile, chain):
    res_dict = get_res_dict(pdbfile, chain)

    cb_lst = []
    tmp_i = 0

    # need to sort to get the sequence correct
    sorted_keys = sorted(res_dict.keys())
    
    for i in sorted_keys:
        if len(res_dict[i]) > 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][-1])
        elif len(res_dict[i]) == 1:
            tmp_i += 1
            cb_lst.append(res_dict[i][0])
            #print tmp_i,i,res_dict[i][0],res_dict[i][-1]
    pdbfile.close()
    return cb_lst

def get_cb_contacts(gapped_cb_lst):
    seqlen = len(gapped_cb_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))

    for i, cb1 in enumerate(gapped_cb_lst):
        if cb1[0] == '-':
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if cb2[0] == '-':
                continue
            diff_vec = cb1 - cb2
            dist_mat[i, j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} <pdb-files> <npz-name> [std in Å]".format(os.path.basename(__file__)))
        sys.exit(1)
    pdb_filename = sys.argv[1]
    npz_name = sys.argv[2]
    if len(sys.argv) > 3:
        std = float(sys.argv[3])
    else:
        std = 1  # One standard deviation is this many Ångström

    STEP = 0.5
    z_per_bin = std/STEP
    z_step = z_per_bin/2
    wanted_bins = (20 - 2)/STEP
    cb_lst = get_cb_coordinates(open(pdb_filename, 'r'), "A")
    contact_mat = get_cb_contacts(cb_lst)
    plen = contact_mat.shape[0]
    dist_mat = np.zeros((plen, plen, 37))

    bins = [i for i in np.arange(2, 2 + (wanted_bins+1)*0.5, 0.5)]
    num_bins = len(bins)

    for i in range(plen):
        for j in range(i, plen):
            if contact_mat[i, j] > 20:
                dist_mat[i, j, 0] = 1
            else:
                dist = contact_mat[i,j]
                ix = np.digitize(dist, bins)
                
                cum_prob = 0
                b_step = 0
                while cum_prob < 0.99:
                    bin_prob = st.norm.cdf(b_step*-z_step)-st.norm.cdf(-z_step*(1+b_step))
                    dist_mat[i, j, np.min([ix+b_step, num_bins-1])] += bin_prob
                    dist_mat[i, j, np.max([ix-b_step, 1])] += bin_prob
                    cum_prob += bin_prob*2
                    b_step += 1

    np.savez_compressed(npz_name, dist=dist_mat)
                    

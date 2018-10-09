import pandas as pd
import numpy as np
import os
import time
import string
import commands
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC
import subprocess
from itertools import chain

flatten = chain.from_iterable

LEFT, RIGHT = 1, -1

def join_ranges(data, offset=0):
    data = sorted(flatten(((start, LEFT), (stop + offset, RIGHT))
            for start, stop in data))
    c = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0:
            yield x, value - offset
            



def get_topidp(seq):
    
    if len(seq) == 0:
        return np.nan
    
    top_idp_avg = 0
        
    for a in top_idp.keys():
        top_idp_avg += (seq.count(a)*top_idp[a])

    top_idp_avg = str(float(top_idp_avg) / float(len(seq)))
    
    #print seq, top_idp_avg
    return top_idp_avg


def get_hessa(seq):
    if len(seq) == 0:
        return np.nan
    
    hessa_avg = 0
        
    for a in hessa.keys():
        hessa_avg += (seq.count(a)*hessa[a])

    hessa_avg = str(float(hessa_avg) / float(len(seq)))
    
    return hessa_avg

def get_ss_scale(seq,dic_scale):
    
    if len(seq) == 0:
        return np.nan
    
    ss_avg = 0
        
    for a in dic_scale.keys():
        ss_avg += (seq.count(a)*dic_scale[a])

    ss_avg = str(float(ss_avg) / float(len(seq)))
    
    return ss_avg

    

    
#~ def parse_iupred(data_file):
    #~ # parse the iupred_multi output
    #~ ps = filter(None, open(data_file).read().split("# Prediction output "))

    #~ iupred_dic = {}
    #~ for p in ps:
        #~ query_id = ""
        #~ values = []

        #~ lines = filter(None,p.split("\n"))
        #~ for line in lines:
            #~ if line[0] == "#":
                #~ query_id = line.split()[1]
            #~ else:
                #~ values += [float(line.split()[-1])]

        #~ # calculate the disorder as the fraction of disordered (>0.5) residues
        #~ if query_id != "":
            #~ '''
            #~ diso_aa = 0
            #~ for a in values:
                #~ if a > 0.5: diso_aa += 1
            
            #~ iupred_dic[query_id] = float(diso_aa)/float(len(values))
            #~ '''
            
            #~ iupred_dic[query_id] = values

    #~ return iupred_dic
            

'''    
def parse_fasta_x(input_file,c_dic):
    # parse the output
    ret_dic = {}
    gene_id = ""
    seq = ""
    x_count = 0
    
    with open(input_file) as output:
        for line in output:
            
            if len(line.strip()) > 0:
                
                if line[0] == ">":
                    
                    gene_id = line[1:].split()[0]
                    seq = ""
                else:
                    seq += line
            else:
                if seq:
                    print gene_id
                    if gene_id in c_dic:
                        final_seq = ""
                        
                        #print c_dic[gene_id]
                        for c_set in c_dic[gene_id]:
                            final_seq += seq[c_set[0]:c_set[1]]
                            ret_dic[gene_id] = float(final_seq.count("x"))/float(len(final_seq))

                            
                            print final_seq
                    seq = ""

        return ret_dic
'''

def parse_fasta_x(input_file,c_dic):

    ret_dic = {}

    recz = SeqIO.parse(input_file, "fasta")
    for r in recz:

        trimmed_id = r.id #.split(".")[0].split("_")[0]

        if trimmed_id in c_dic:
            final_seq = ""
            for c_set in c_dic[trimmed_id]:
                final_seq += r.seq[c_set[0]:c_set[1]]

            result = final_seq.count("x")    
            ret_dic[r.id] = result
            
            #print r.id, result
            #print final_seq

    return ret_dic

def aa_freq(seq,aa):
    if len(seq) == 0:
        return np.nan
    #return float(seq.count(aa)) / float(len(seq))
    return seq.count(aa)

      
def aa_count(seq,aa):
    if len(seq) == 0:
        return np.nan
    return seq.count(aa)



def get_regions_coords(recs,ids):
    
    dic_results = {}
    dic_results["full"] = {}
    dic_results["domains_shared"] = {}
    dic_results["domains_other"] = {}
    dic_results["linkers"] = {}
    
    for r in ids:
        uniprot_id = recs[r].id.split(".")[0].split("_")[0]
        
        coords_all = [(0,len(recs[r]))]
        #print uniprot_id

        ## EXTRACT LINKERS COORDINATES
        df_coords_domains_selected = df_regions_sel.loc[df_regions_sel.uniprot_acc == uniprot_id]

        ranges = []
        for index,row in df_coords_domains_selected[['seq_start','seq_end']].iterrows():
            ranges += [[row.seq_start,row.seq_end]]

        ranges = list(join_ranges(ranges))
        subl = [item for sublist in ranges for item in sublist]


        coords_linkers = []

        if subl[0] != 0:
            subl = [0] + subl + [len(recs[r])]
            for p in range(len(subl)):
                if p % 2 == 0:
                    coords_linkers += [(subl[p],subl[p+1])]
        else:
            subl = subl + [len(recs[r])]
            for p in range(len(subl)-1):
                if p % 2 == 1:
                    coords_linkers += [(subl[p],subl[p+1])] 



        ## EXTRACT SHARED DOMAINS COORDINATES
        df_coords_domains_selected_shared = df_coords_domains_selected.loc[df_coords_domains_selected.is_shared == True]

        ranges = []
        for index,row in df_coords_domains_selected_shared[['seq_start','seq_end']].iterrows():
            ranges += [[row.seq_start,row.seq_end]]

        ranges = list(join_ranges(ranges))
        subl = [item for sublist in ranges for item in sublist]

        coords_shared_domains = []
        for p in range(len(subl)):
            if p % 2 == 0:
                coords_shared_domains += [(subl[p],subl[p+1])]


        ## EXTRACT OTHER DOMAINS COORDINATES
        df_coords_domains_selected_other = df_coords_domains_selected.loc[df_coords_domains_selected.is_shared == False]

        ranges = []
        for index,row in df_coords_domains_selected_other[['seq_start','seq_end']].iterrows():
            ranges += [[row.seq_start,row.seq_end]]

        ranges = list(join_ranges(ranges))
        subl = [item for sublist in ranges for item in sublist]

        coords_other_domains = []
        for p in range(len(subl)):
            if p % 2 == 0:
                coords_other_domains += [(subl[p],subl[p+1])]


        #dic_results[r] = coords_all, coords_shared_domains, coords_other_domains, coords_linkers
        dic_results["full"][r] = coords_all
        dic_results["domains_shared"][r] = coords_shared_domains
        dic_results["domains_other"][r] = coords_other_domains
        dic_results["linkers"][r] = coords_linkers
    
    return dic_results
    
    
    



############################################################################
##    MAIN
############################################################################

dir='/scratch2/arne/proks-vs-euk/'

if  (not os.path.isdir(dir)):
    dir='/pfs/nobackup/home/w/wbasile/proks_euks/'
        
data_dir = dir + "/data/"
input_dir = dir + "/results/uniprot_sequences_groups/"
annotations_dir = dir + "/results/uniprot_pfam_annotations/"

# MAIN

# Top-IDP scale
top_idp = {'A':0.06, 'C' :  0.02, 'D' : 0.192, 'E' : 0.736,
    'F' :  -0.697, 'G' : 0.166, 'H':0.303, 'I' :  -0.486,
    'K' : 0.586, 'L' :  -0.326, 'M': -0.397, 'N' : 0.007,
    'P' : 0.987, 'Q' : 0.318, 'R' : 0.180, 'S':  0.341,
    'T' : 0.059, 'V' :  -0.121, 'W':  -0.884, 'Y' : -0.510}

# Hessa scale
hessa = {
'C':-0.13,
'D':3.49,
'S':0.84,
'Q':2.36,
'K':2.71,
'W':0.3,
'P':2.23,
'T':0.52,
'I':-0.6,
'A':0.11,
'F':-0.32,
'G':0.74,
'H':2.06,
'L':-0.55,
'R':2.58,
'M':-0.1,
'E':2.68,
'N':2.05,
'Y':0.68,
'V':-0.31
}

dic_aa = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS',
              'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 
              'V': 'VAL', 'Y': 'TYR'}




df_scales = pd.read_csv(dir + "/data/scales_and_slopes.tsv",sep="\t")
df_scales.AA = df_scales.AA.apply(string.upper)

dic_aa_inv = {}
for a in dic_aa:
    dic_aa_inv[dic_aa[a]] = a
    
df_scales.AA = df_scales.AA.map(dic_aa_inv)

dic_ss_alpha = df_scales.set_index("AA").to_dict()["Alpha"]
dic_ss_beta = df_scales.set_index("AA").to_dict()["Beta"]
dic_ss_coil = df_scales.set_index("AA").to_dict()["Coil"]
dic_ss_turn = df_scales.set_index("AA").to_dict()["Turn"]



aas = ['A',
 'C',
 'D',
 'E',
 'F',
 'G',
 'H',
 'I',
 'K',
 'L',
 'M',
 'N',
 'P',
 'Q',
 'R',
 'S',
 'T',
 'V',
 'W',
 'Y']



#~ input_dir = dir + "/data/untitled_folder/"
#~ annotations_dir = dir + "/data/test_res/"

all_sets = ["full", "domains_shared", "domains_other", "linkers"]



if not os.path.exists(annotations_dir):
    os.makedirs(annotations_dir)

    
# load the list of PFAM domains that are shared by at least 5 bacteria and 5 euks
out_domain_ids_filename = dir+"bin/pfam_ids_orthologs_10.list"
shared_domains_pfam_ids = set(filter(None, open(out_domain_ids_filename).read().split("\n")))

pfam_regions_selected_filename = data_dir + "/pfam/Pfam-A.regions.uniprot.selected_proteins.selected_fields.feather"
print "loading coords file " + pfam_regions_selected_filename  + " ...",
df_regions = pd.read_feather(pfam_regions_selected_filename)#[['uniprot_acc', u'pfamA_acc', u'seq_start', u'seq_end']]
#df_regions = pd.read_csv(dir + "/data/pfam/PfamA.regions.uniprot.test.tsv", sep="\t")#[['uniprot_acc', u'pfamA_acc', u'seq_start', u'seq_end']]

print "done"

    
file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        file_list += [input_dir + f]
        



for f in file_list:
    
    f_name = f.split("/")[-1]
    out_annotation_file = annotations_dir + f_name +"_annotation_"
    
    if not os.path.exists(out_annotation_file):
        
        # check the existence of these three data files (produced separately)
        iupred_long_data_file = f +  ".data_iupred_long"
        iupred_short_data_file = f +  ".data_iupred_short"
        seg_data_file = f +  ".data_seg"
        
        proceed = True
        
        if not os.path.exists(iupred_long_data_file):
            proceed = False
 
        if not os.path.exists(iupred_short_data_file):
            proceed = False
 
        if not os.path.exists(seg_data_file):
            proceed = False
 
        
        if proceed == True:
            
            print "processing " + f + " to " + out_annotation_file
        
            # create a placeholder file to prevent other processes to work on this file
            open(out_annotation_file, 'a').close()
            
            recs = SeqIO.index(f,"fasta")
            
            # for each sequence, obtain 4 sets of coordinates (full, domains_shared, domains_other, linkers)
            valid_ids = []
            counter_dic = {}
            for r in recs:
                rid = recs[r].id
                counter_id = rid.split(".")[0].split("_")[0]
                valid_ids += [counter_id]
                counter_dic[counter_id] = rid

            print "filtering coords file ...",
            df_regions_sel = df_regions.loc[df_regions.uniprot_acc.isin(valid_ids)]#[['uniprot_acc', u'pfamA_acc', u'seq_start', u'seq_end']]
            print "done"
            
            # this is the set of uniprot ids present in both the coordinates file ( == orthologs 5) and in uniprot file
            final_ids = set([counter_dic[r] for r in df_regions_sel.uniprot_acc])
            
            df_regions_sel["is_shared"] = df_regions_sel.pfamA_acc.isin(shared_domains_pfam_ids)
            
            print "parsing coords file ...",
            coords_dic = get_regions_coords(recs,final_ids)
            print "done"
                 
            
            # create 4 dataframes (full, domains_shared, domains_other, linkers)
            dataframes = {}
            for ty in all_sets:
                dataframes[ty] = pd.DataFrame()
                dataframes[ty]["query_id"] = list(final_ids)
                
                
            # prepare dictionaries to hold the different annotations
            dic_sequence = {}
            for ty in all_sets:
                dic_sequence[ty] = {}
                
            
            # annotate each protein
            for rid in final_ids:
                rec = recs[rid]
                rec_seq = str(rec.seq)
                
                # extract the sequences
                for ty in all_sets:    
                    dic_sequence[ty][rid] = ""
                    for c_set in coords_dic[ty][rid]:
                        dic_sequence[ty][rid] += rec_seq[c_set[0]:c_set[1]]
                    
            
            for ty in all_sets:

                dataframes[ty]["seq"] = dataframes[ty]["query_id"].map(dic_sequence[ty])
            
                # discard all records without sequence (for ex. proteins with no "other" domains)
                dataframes[ty] = dataframes[ty].dropna(subset=["seq"])

                # annotate all properties
                dataframes[ty]["length_translation"] = dataframes[ty]["seq"].apply(len)
                dataframes[ty]["top-idp"] = dataframes[ty]["seq"].apply(get_topidp)
                dataframes[ty]["hessa"] = dataframes[ty]["seq"].apply(get_hessa)
                
                # AA frequencies
                for aa in aas:
                    dataframes[ty]["freq_" + aa] = dataframes[ty]["seq"].apply(aa_freq, args = (aa,)) 

                # AA counts
                #for aa in aas:
                #    df["count_" + aa] = df.seq_translation.apply(aa_count, args = (aa,))

                ## add SS scales
                for ss_type, d_scale in zip(["alpha", "beta", "coil", "turn"],[dic_ss_alpha, dic_ss_beta, dic_ss_coil, dic_ss_turn]):
                    dataframes[ty]["ss_" + ss_type] = dataframes[ty]["seq"].apply(get_ss_scale, args = (d_scale,))
                    
                
                                    
                # add iupred
                #print "iupred long"
                dic_iupred_long = parse_fasta_x(iupred_long_data_file,coords_dic[ty])
                dataframes[ty]['iupred_long'] = dataframes[ty]['query_id'].map(dic_iupred_long)
                
                #print "iupred short"
                dic_iupred_short = parse_fasta_x(iupred_short_data_file,coords_dic[ty])
                dataframes[ty]['iupred_short'] = dataframes[ty]['query_id'].map(dic_iupred_short)

                
                # add SEG
                #print "SEG"
                dic_seg = parse_fasta_x(seg_data_file,coords_dic[ty])
                dataframes[ty]['seg'] = dataframes[ty]['query_id'].map(dic_seg)
                
                
                # export
                columns = ["query_id",  "length_translation", "top-idp", "hessa", "iupred_long", "iupred_short", "seg"]
                columns += ["freq_" + aa for aa in aas]
                columns += ["ss_alpha", "ss_beta", "ss_coil", "ss_turn"]
                
                dataframes[ty][columns].to_csv(out_annotation_file + ty + ".csv", index=False, header = False)
                
                # export the header separately
                with open(annotations_dir + f_name +"_annotation.header", "w") as houtf:
                    houtf.write(",".join(columns) + "\n")
                    
        
                #print 
                



# coding: utf-8

# In[1]:

# create a dataframe containing only the queries
# each row is a query orf, with its corresponding protein, translation and all the intrinsic properties 
# (length, gc, disorder, TM, hydrophobicity (hessa), top-idp, mreps, SEG, volatility...)

#from settings import *
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC
import subprocess


################################
ty = "domains_exclusive_euk"
###############################


data_dir = "../data/"
input_dir = "../results/pfam-"+ty+"/"
annotations_dir = "../results/pfam-"+ty+"_annotations/"

if not os.path.exists(annotations_dir):
    os.makedirs(annotations_dir)


# In[2]:




# In[3]:

# Top-IDP scale
top_idp = {'A':0.06, 'C' :  0.02, 'D' : 0.192, 'E' : 0.736,
    'F' :  -0.697, 'G' : 0.166, 'H':0.303, 'I' :  -0.486,
    'K' : 0.586, 'L' :  -0.326, 'M': -0.397, 'N' : 0.007,
    'P' : 0.987, 'Q' : 0.318, 'R' : 0.180, 'S':  0.341,
    'T' : 0.059, 'V' :  -0.121, 'W':  -0.884, 'Y' : -0.510}

def get_topidp(seq):
    
    if len(seq) == 0:
        return np.nan
    
    top_idp_avg = 0
        
    for a in top_idp.keys():
        top_idp_avg += (seq.count(a)*top_idp[a])

    top_idp_avg = str(float(top_idp_avg) / float(len(seq)))
    
    return top_idp_avg


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

def get_hessa(seq):
    if len(seq) == 0:
        return np.nan
    
    hessa_avg = 0
        
    for a in hessa.keys():
        hessa_avg += (seq.count(a)*hessa[a])

    hessa_avg = str(float(hessa_avg) / float(len(seq)))
    
    return hessa_avg



dic_aa = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS',
              'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 
              'V': 'VAL', 'Y': 'TYR'}


import string

df_scales = pd.read_csv("../data/scales_and_slopes.tsv",sep="\t")
df_scales.AA = df_scales.AA.apply(string.upper)

dic_aa_inv = {}
for a in dic_aa:
    dic_aa_inv[dic_aa[a]] = a
    
df_scales.AA = df_scales.AA.map(dic_aa_inv)

dic_ss_alpha = df_scales.set_index("AA").to_dict()["Alpha"]
dic_ss_beta = df_scales.set_index("AA").to_dict()["Beta"]
dic_ss_coil = df_scales.set_index("AA").to_dict()["Coil"]
dic_ss_turn = df_scales.set_index("AA").to_dict()["Turn"]


def get_ss_scale(seq,dic_scale):
    
    if len(seq) == 0:
        return np.nan
    
    ss_avg = 0
        
    for a in dic_scale.keys():
        ss_avg += (seq.count(a)*dic_scale[a])

    ss_avg = str(float(ss_avg) / float(len(seq)))
    
    return ss_avg

    

    
def parse_iupred(data_file):
    # parse the iupred_multi output
    ps = filter(None, open(data_file).read().split("# Prediction output "))

    iupred_dic = {}
    for p in ps:
        query_id = ""
        values = []

        lines = filter(None,p.split("\n"))
        for line in lines:
            if line[0] == "#":
                query_id = line.split()[1]
            else:
                values += [float(line.split()[-1])]

        # calculate the disorder as the fraction of disordered (>0.5) residues
        if query_id != "":
            diso_aa = 0
            for a in values:
                if a > 0.5: diso_aa += 1

            iupred_dic[query_id] = float(diso_aa)/float(len(values))

    return iupred_dic
            
    
    
def do_seg(input_file):
    output = subprocess.check_output('./seg '+ input_file + " -x", shell=True,stderr=subprocess.PIPE)

    # parse the output
    seg_dic = {}
    gene_id = ""
    seq = ""
    x_count = 0
    for line in output.split("\n"):
        if len(line) > 0:
            if line[0] == ">":
                gene_id = line[1:].split()[0]
                seq = ""
            else:
                seq += line
        else:
            if seq:
                seg_dic[gene_id] = float(seq.count("x"))/float(len(seq))
                seq = ""

    return seg_dic


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


def aa_freq(seq,aa):
    if len(seq) == 0:
        return np.nan
    return float(seq.count(aa)) / float(len(seq))

      
def aa_count(seq,aa):
    if len(seq) == 0:
        return np.nan
    return seq.count(aa)
      


def annotate_genome(f):

    if not os.path.exists("./stop"):

        print f

        f_name = f.split("/")[-1]

        out_annotation_file = annotations_dir + f_name +"_annotation.csv"
        

        if not os.path.exists(out_annotation_file):
            
            # create an ampty file to prevent the other processes to work on the same genome
            cmd = "touch " + out_annotation_file
            os.system(cmd)

            protein_recs = SeqIO.to_dict(SeqIO.parse(f,"fasta"))

            recs = []
            for k in protein_recs:
                rec = {}
                rec["query_id"] = k
                rec["seq_translation"] = str(protein_recs[k].seq)
                recs += [rec]

            df = pd.DataFrame(recs)

            # annotate all properties
            df["length_translation"] = df.seq_translation.apply(len)
            df["top-idp"] = df["seq_translation"].apply(get_topidp)
            df["hessa"] = df["seq_translation"].apply(get_hessa)

            # AA frequencies
            for aa in aas:
                df["freq_" + aa] = df.seq_translation.apply(aa_freq, args = (aa,)) 

            # AA counts
            for aa in aas:
                df["count_" + aa] = df.seq_translation.apply(aa_count, args = (aa,))

            ## add SS scales
            for ss_type, d_scale in zip(["alpha", "beta", "coil", "turn"],[dic_ss_alpha, dic_ss_beta, dic_ss_coil, dic_ss_turn]):
                df["ss_" + ss_type] = df.seq_translation.apply(get_ss_scale, args = (d_scale,))

                
            # disorder
            print str(f),"Computing IUpred long"
            iupred_long_data_raw_file = f +  "data_iupred_long" 

            if os.path.exists(iupred_long_data_raw_file):
                st = os.stat(iupred_long_data_raw_file)
                if st.st_size == 0:
                    os.system("rm " + iupred_long_data_raw_file)

            if not os.path.exists(iupred_long_data_raw_file):
                cmd = "./iupred_multi " + f + " long > " + iupred_long_data_raw_file
                print cmd
                os.system(cmd)
            iupred_dic = parse_iupred(iupred_long_data_raw_file)
            df['iupred_long'] = df['query_id'].map(iupred_dic)
            #os.system("rm " + iupred_long_data_raw_file)


            print str(f),"Computing IUpred short"
            iupred_short_data_raw_file = f +  "data_iupred_short"

            if os.path.exists(iupred_short_data_raw_file):
                st = os.stat(iupred_short_data_raw_file)
                if st.st_size == 0:
                    os.system("rm " + iupred_short_data_raw_file)
            
            if not os.path.exists(iupred_short_data_raw_file):
                cmd = "./iupred_multi " + f + " short > " + iupred_short_data_raw_file
                print cmd
                os.system(cmd)
            iupred_dic = parse_iupred(iupred_short_data_raw_file)
            df['iupred_short'] = df['query_id'].map(iupred_dic)
            #os.system("rm " + iupred_short_data_raw_file)

            #df['iupred_long'] = np.nan
            #df['iupred_short'] = np.nan


            # SEG
            print str(f),"Computing SEG"
            seg_dic = do_seg(f)
            df["seg"] = df["query_id"].map(seg_dic)


            # export
            columns = ["query_id",  "length_translation", "top-idp", "hessa", "iupred_long", "iupred_short", "seg"]
            columns += ["freq_" + aa for aa in aas]
            columns += ["count_" + aa for aa in aas]
            columns += ["ss_alpha", "ss_beta", "ss_coil", "ss_turn"]
            df[columns].dropna(subset=["hessa"]).to_csv(out_annotation_file, index = False)
            #~ df.to_csv(out_annotation_file, index = False)

        
#pool = Pool()
#res = pool.map(annotate_genome,file_list)

file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        file_list += [input_dir + f]
        
#~ file_list = ["../data/pfam/pfam_domains_exclusive_bac.fsa"]

for f in file_list:
    annotate_genome(f)




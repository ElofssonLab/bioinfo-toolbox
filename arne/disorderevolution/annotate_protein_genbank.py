#!/usr/bin/env python3
# coding: utf-8
import pandas as pd
import numpy as np
import os
import re
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC
import subprocess
from Bio import SwissProt
#from Bio.SwissProt import KeyWList



def get_topidp(seq):
    
    if len(seq) == 0:
        return np.nan
    
    top_idp_avg = 0
        
    for a in top_idp.keys():
        top_idp_avg += (seq.count(a)*top_idp[a])

    top_idp_avg = str(float(top_idp_avg) / float(len(seq)))
    
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
    segexec=dir+"/bin/seg  "
    seg_dic = {}
    gene_id = ""
    seq = ""
    x_count = 0
    #print(segexec + input_file + " -x")
    try:
        output = subprocess.check_output(segexec + input_file + " -x", shell=True,stderr=subprocess.PIPE)
        #print (output)
    except:
        with open(input_file, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                ID=re.sub(r'.*\|','',record.id)
                seg_dic[ID]=0
                # parse the output
    for line in output.split("\n"):
        if len(line) > 0:
            if line[0] == ">":
                gene_id = line[1:].split()[0]
                seq = ""
            else:
                seq += line
        else:
            if seq:
                ID=re.sub(r'.*\|','',gene_id)
                seg_dic[ID] = float(seq.count("x"))/float(len(seq)+1.e-20)
                seq = ""
    return seg_dic

def aa_freq(seq,aa):
    if len(seq) == 0:
        return np.nan
    return float(seq.count(aa)) / float(len(seq))

      
def aa_count(seq,aa):
    if len(seq) == 0:
        return np.nan
    return seq.count(aa)
      

def parse_fasta_x(input_file):
    ret_dic = {}
    recz = SeqIO.parse(input_file, "fasta")
    for r in recz:
        result = float(r.seq.count("x")) / float(len(r))
        id=re.sub(r'.*\|','',r.id)
        ret_dic[id] = result
        #print (r,result)
    return ret_dic


def parse_uniprot(input_file):
    dic_pfam = {}
    dic_dom = {}
    dic_king = {}
    # probably faster/easier to use the XML parser directly
    #print (input_file)
    handle = open(input_file)
    for record in SwissProt.parse(handle):
        #print (record)
        #print (record.entry_name)
        #print (record.cross_references)
        entry =record.entry_name
        id=entry
        dic_pfam[id]=''
        dic_dom[id]=0
        dic_king[id]='Unique'
        for  db in record.cross_references:
            if ( db[0] == "Pfam"):
                dic_pfam[id]=dic_pfam[id]+db[1]+";"
                dic_dom[id]+=1
                if (db[1] in shared_domains.keys()):
                    dic_king[id]="Shared"
        if (dic_dom[id]==0):
            dic_king[id]='None'
    return dic_pfam,dic_dom,dic_king



def annotate_genome(f):

    if not os.path.exists("./stop"):
        # check the existence of these three data files (produced separately)
        iupred_long_data_file = f +  ".data_iupred_long"
        iupred_short_data_file = f +  ".data_iupred_short"
        iupred04_long_data_file = f +  ".data_iupred_0.4_long"
        iupred04_short_data_file = f +  ".data_iupred_0.4_short"
        uniprot_data_file = re.sub(r'\.fasta','.txt',f)

        proceed = True
        
        if not os.path.exists(iupred_long_data_file):
            print ("Missing IUPRED-long")
            proceed = False
 
        if not os.path.exists(iupred_short_data_file):
            print ("Missing IUPRED-short")
            proceed = False

        if not os.path.exists(iupred04_long_data_file):
            print ("Missing IUPRED04-long")
            proceed = False
 
        if not os.path.exists(iupred04_short_data_file):
            print ("Missing IUPRED04-short")
            proceed = False

        if proceed == True:
            #print ("Trying:",f)

            f_name = f.split("/")[-1]

            out_annotation_file = output_dir + f_name +"_annotation.csv"
            
            if not os.path.exists(out_annotation_file):
                
                # create an ampty file to prevent the other processes to work on the same genome
                cmd = "touch " + out_annotation_file
                os.system(cmd)
                protein_recs = SeqIO.to_dict(SeqIO.parse(f,"fasta"))
                recs = []
                for k in protein_recs:
                    id=re.sub(r'.*\|','',k)
                    rec = {}
                    rec["query_id"] = id
                    rec["seq"] = str(protein_recs[k].seq)
                    recs += [rec]

                df = pd.DataFrame(recs)

                # annotate all properties
                df["length"] = df.seq.apply(len)
                df["top-idp"] = df["seq"].apply(get_topidp)
                df["hessa"] = df["seq"].apply(get_hessa)


                # AA frequencies
                for aa in aas:
                    df["freq_" + aa] = df.seq.apply(aa_freq, args = (aa,)) 

                # AA counts
                #for aa in aas:
                #    df["count_" + aa] = df.seq.apply(aa_count, args = (aa,))

                ## add SS scales
                for ss_type, d_scale in zip(["alpha", "beta", "coil", "turn"],[dic_ss_alpha, dic_ss_beta, dic_ss_coil, dic_ss_turn]):
                    df["ss_" + ss_type] = df.seq.apply(get_ss_scale, args = (d_scale,))

                    
                # disorder

                # add iupred
                print ("iupred long")
                dic_iupred_long = parse_fasta_x(iupred_long_data_file)
                df['iupred_long'] = df['query_id'].map(dic_iupred_long)
                #print (dic_iupred_long,df['iupred_long'],df['query_id'])
                
                print ("iupred short")
                dic_iupred_short = parse_fasta_x(iupred_short_data_file)
                df['iupred_short'] = df['query_id'].map(dic_iupred_short)

                print ("iupred long 0.4")
                dic_iupred04_long = parse_fasta_x(iupred04_long_data_file)
                df['iupred04_long'] = df['query_id'].map(dic_iupred04_long)
                #print (dic_iupred04_long,df['iupred04_long'],df['query_id'])
                
                print ("iupred short 0.4")
                dic_iupred04_short = parse_fasta_x(iupred04_short_data_file)
                df['iupred04_short'] = df['query_id'].map(dic_iupred04_short)


                # SEG
                print (str(f),"Computing SEG")
                seg_dic = do_seg(f)
                df["seg"] = df["query_id"].map(seg_dic)


                # Parsing uniprot - today only Pfam
                #print "Parse uniprot"
                (dic_pfam,dic_numdoms,dic_kingdom) = parse_uniprot(uniprot_data_file)
                df["Pfam"] = df['query_id'].map(dic_pfam)
                df["NumDoms"] = df['query_id'].map(dic_numdoms)
                df["PfamType"] = df['query_id'].map(dic_kingdom)

                # export
                columns = ["query_id",  "length", "top-idp", "iupred_long", "iupred_short","iupred04_long", "iupred04_short","seg","ss_alpha", "ss_beta", "ss_coil", "ss_turn","hessa","Pfam","NumDoms","PfamType"]
                columns += ["freq_" + aa for aa in aas]
                df[columns].to_csv(out_annotation_file, index = False)

if __name__ == '__main__':

	# This is where we start
	
	output_dir=':~/git/paper_proks_vs_euks_proteins/data/'
        input_dir = '/scratch2/arne/genbank/'
	
	if not os.path.exists(output_dir):
	    os.makedirs(output_dir)
	
	
	# Top-IDP scale
	top_idp = {'A':0.06, 'C' :  0.02, 'D' : 0.192, 'E' : 0.736,
	    'F' :  -0.697, 'G' : 0.166, 'H':0.303, 'I' :  -0.486,
	    'K' : 0.586, 'L' :  -0.326, 'M': -0.397, 'N' : 0.007,
	    'P' : 0.987, 'Q' : 0.318, 'R' : 0.180, 'S':  0.341,
	    'T' : 0.059, 'V' :  -0.121, 'W':  -0.884, 'Y' : -0.510}
	
	
	
	# Hessa scale
	hessa = { 'C':-0.13, 'D':3.49, 'S':0.84, 'Q':2.36, 'K':2.71,
	          'W':0.3, 'P':2.23, 'T':0.52, 'I':-0.6, 'A':0.11, 'F':-0.32,
	          'G':0.74, 'H':2.06, 'L':-0.55, 'R':2.58, 'M':-0.1, 'E':2.68,
	          'N':2.05, 'Y':0.68, 'V':-0.31 }
	
	
	
	
	dic_aa = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY',
	          'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M':
	          'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO',
	          'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V':
	          'VAL', 'Y': 'TYR'}
	
	
	import string
	
	df_scales = pd.read_csv(data_dir + "scales_and_slopes.tsv",sep="\t")
	df_scales.AA = df_scales.AA.apply(string.upper)
	
	dic_aa_inv = {}
	for a in dic_aa:
	    dic_aa_inv[dic_aa[a]] = a
	    
	df_scales.AA = df_scales.AA.map(dic_aa_inv)
	
	dic_ss_alpha = df_scales.set_index("AA").to_dict()["Alpha"]
	dic_ss_beta = df_scales.set_index("AA").to_dict()["Beta"]
	dic_ss_coil = df_scales.set_index("AA").to_dict()["Coil"]
	dic_ss_turn = df_scales.set_index("AA").to_dict()["Turn"]
	
	
	    
	
	
	aas = ['A','C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
	       'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	
	nucl = ['A','C','T','G']

	GC = ['C','G']
	            
	        
	
	
	file_list = []
	for f in os.listdir(input_dir):
	    if f.endswith(".gbk.gz"):
	        file_list += [input_dir + f]
	
	#annotate_genome('data/proteomes/UP000001554_7739.fasta')
        #sys.exit()
	for f in file_list:
	    print (f)
	    try:
	        annotate_genome(f)
	    except:
	        print ("ERROR on " + f)





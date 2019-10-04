
# coding: utf-8

import pandas as pd
import numpy as np
import os
import re
import sys
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


dir='/scratch2/arne/annotate_uniprot_proteomes/'

if  (not os.path.isdir(dir)):
    dir='/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/'

data_dir = dir+"/data/"
input_dir = data_dir + "proteomes/"

output_dir = dir+"/results/extended/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)


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


global nucleotides,codons
nucleotides= ["A","C","T","G"]
codons=[]
nucleotidepos=[]
for one in nucleotides:
    for pos in ["1","2","3"]:
        nucleotidepos += [one+pos]
        for two in nucleotides:
            for three in nucleotides:
                codons+=[one+two+three]

                
# Hessa scale


hessa = {'C':-0.13, 'D':3.49, 'S':0.84, 'Q':2.36,
'K':2.71, 'W':0.3, 'P':2.23, 'T':0.52, 'I':-0.6, 'A':0.11, 'F':-0.32,
'G':0.74, 'H':2.06, 'L':-0.55, 'R':2.58, 'M':-0.1, 'E':2.68, 'N':2.05,
'Y':0.68, 'V':-0.31 }

def get_hessa(seq):
    if len(seq) == 0:
        return np.nan
    hessa_avg = 0
    for a in hessa.keys():
        hessa_avg += (seq.count(a)*hessa[a])
    hessa_avg = str(float(hessa_avg) / float(len(seq)))
    return hessa_avg



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

def parse_scampi(data_file):
    scampi_i = {}
    scampi_m = {}
    scampi_o = {}
    rec = SeqIO.parse(data_file, "fasta")
    for r in rec:
        i = float(r.seq.count("I")) / float(len(r))
        o = float(r.seq.count("O")) / float(len(r))
        m = float(r.seq.count("M")) / float(len(r))
        id=re.sub(r'.*\|','',r.id)
        scampi_i[id] = i
        scampi_o[id] = o
        scampi_m[id] = m

    return scampi_i,scampi_m,scampi_o
            
    
def parse_scampi_seq(scampi_file,fasta_file):
    # This is not ready..

    scampi = SeqIO.to_dict(SeqIO.parse(scampi_file, "fasta"))
    seq = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    freq={}
    sumtype={}
    for type in ["I","O","M"]:
        for aa in  ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
            freq[type+"-"+aa]={}

    for r in scampi:
        
        for type in ["I","O","M"]:
            sumtype[type]=0
        if (scampi[r].seq.count("M")>0):
            for type in ["I","O","M"]:
                sumtype[type]=scampi[r].seq.count(type) 
            print (seq[r].seq,len(seq[r]),sumtype)
            print (scampi[r].seq,len(scampi[r]))
            for i in range(0,len(seq[r])):
                type=scampi[r].seq[i]
                aa=seq[r].seq[i]
                freq[type+"-"+aa]+=1/sumtype[type]

            
    return freq
            
    
    
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


aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def aa_freq(seq,aa):
    if len(seq) == 0:
        return np.nan
    return float(seq.count(aa)) / float(len(seq))

def gc_freq(seq,pos):
    if len(seq) == 0:
        return np.nan
    tempseq=''
    for i in range(pos-1,len(seq),3):
        tempseq+=seq[i]
    return float(tempseq.count('G')+tempseq.count('C')) / float(len(tempseq))

def nucl_freq(seq,pos,nucl):
    if len(seq) == 0:
        return np.nan
    nucfreq=0.
    tempseq=''
    for i in range(pos-1,len(seq),3):
        tempseq+=seq[i]
    nucfreq=float(tempseq.count(nucl)) / float(len(tempseq))
    return nucfreq

def codon_freq(seq,codon):
    if len(seq) == 0:
        return np.nan
    codonfreq=0.
    for i in range(0,len(seq),3):
        tempcodon=seq[i:i+3]
        if (tempcodon == codon):
            codonfreq+=3./float(len(seq))
    return codonfreq


      
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
        scampi_data_file = f +  ".scampi"
        uniprot_data_file = re.sub(r'\.fasta','.txt',f)
        DNA_fasta_file = re.sub(r'\.fasta','_DNA.fasta',f)

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

        if not os.path.exists(scampi_data_file):
            print ("Missing SCAMPI")
            proceed = False
            
        if proceed == True:
            print ("Trying:",f)

            f_name = f.split("/")[-1]

            out_annotation_file = output_dir + f_name +"_annotation.csv"
            
            if not os.path.exists(out_annotation_file):
                
                # create an ampty file to prevent the other processes to work on the same genome
                cmd = "touch " + out_annotation_file
                os.system(cmd)
                protein_recs = SeqIO.to_dict(SeqIO.parse(f,"fasta"))
                recs = []
                for k in protein_recs:
                    longid=re.sub(r'.*\|','',k)
                    #id=re.sub(r'\>[a-z]+\|','',re.sub(r'\|.*','',k))
                    id=re.split('\|',k)
                    #print (id,longid,k)
                    rec = {}
                    rec["query_id"] = longid
                    rec["shortid"] = id[1]
                    rec["seq"] = str(protein_recs[k].seq)
                    recs += [rec]

                df_prot = pd.DataFrame(recs)

                # annotate all properties
                df_prot["length"] = df_prot.seq.apply(len)
                df_prot["top-idp"] = df_prot["seq"].apply(get_topidp)
                df_prot["hessa"] = df_prot["seq"].apply(get_hessa)


                # AA frequencies
                for aa in aas:
                    df_prot["freq_" + aa] = df_prot.seq.apply(aa_freq, args = (aa,)) 

                for ss_type, d_scale in zip(["alpha", "beta", "coil", "turn"],[dic_ss_alpha, dic_ss_beta, dic_ss_coil, dic_ss_turn]):
                    df_prot["ss_" + ss_type] = df_prot.seq.apply(get_ss_scale, args = (d_scale,))

                    # AA counts
                #for aa in aas:
                #    df["count_" + aa] = df_prot.seq.apply(aa_count, args = (aa,))

                ## add SS scales

                    
                # disorder

                # add iupred
                #print "iupred long"
                dic_iupred_long = parse_fasta_x(iupred_long_data_file)
                df_prot['iupred_long'] = df_prot['query_id'].map(dic_iupred_long)
                #print (dic_iupred_long,df_prot['iupred_long'],df_prot['query_id'])
                
                #print "iupred short"
                dic_iupred_short = parse_fasta_x(iupred_short_data_file)
                df_prot['iupred_short'] = df_prot['query_id'].map(dic_iupred_short)

                #print "iupred long 0.4"
                dic_iupred04_long = parse_fasta_x(iupred04_long_data_file)
                df_prot['iupred04_long'] = df_prot['query_id'].map(dic_iupred04_long)
                #print (dic_iupred04_long,df_prot['iupred04_long'],df_prot['query_id'])
                
                #print "iupred short 0.4"
                dic_iupred04_short = parse_fasta_x(iupred04_short_data_file)
                df_prot['iupred04_short'] = df_prot['query_id'].map(dic_iupred04_short)

                #print "iupred long"
                dic_scampi_i,dic_scampi_m,dic_scampi_o = parse_scampi(scampi_data_file)
                df_prot['scampi_i'] = df_prot['query_id'].map(dic_scampi_i)
                df_prot['scampi_m'] = df_prot['query_id'].map(dic_scampi_m)
                df_prot['scampi_o'] = df_prot['query_id'].map(dic_scampi_o)
                #print (dic_iupred_long,df_prot['iupred_long'],df_prot['query_id'])

                #scampi_freq = parse_scampi_seq(scampi_data_file,f)
                #print (scampi_freq)
                #sys.exit()
                # SEG
                # print str(f),"Computing SEG"
                seg_dic = do_seg(f)
                df_prot["seg"] = df_prot["query_id"].map(seg_dic)



                # GC counts
                dna_recs = SeqIO.to_dict(SeqIO.parse(DNA_fasta_file,"fasta"))
                dnarecs = []
                for k in dna_recs:
                    #id=re.sub(r'\|.*','',re.sub(r'\>[a-z]+\|','',k))
                    id=re.split('\|',k)
                    dnarec = {}
                    #print (k,id[1])
                    dnarec["shortid"] = id[1]
                    dnarec["dnaseq"] = str(dna_recs[k].seq)
                    dnarecs += [dnarec]

                #print ("RECS",len(recs))
                #print ("DNA",len(dnarecs))
                df_dna=pd.DataFrame(dnarecs)

                    #print(df)
                #print(tempdf)
                sum=0
                for i in range(1,4):
                    df_dna["GC"+str(i)] = df_dna.dnaseq.apply(gc_freq,args = (i,))
                    sum+=df_dna["GC"+str(i)]
                df_dna["GCcoding"]= sum/3.

                # Now also indiviual codons
                #nucl=nucl_freq(dnaseq)
                #print (nucl)
                for pos in range(1,4):
                    for nucl in nucleotides:
                        df_dna[nucl+str(pos)]=df_dna.dnaseq.apply(nucl_freq, args = (pos,nucl,))
                for nucl in nucleotides:
                    df_dna[nucl]=(df_dna[nucl+"1"]+df_dna[nucl+"2"]+df_dna[nucl+"3"])/3
                #print (df_dna)
                #sys.exit()
                # And codons
                #codonfreq={}
                #codonfreq=codon_freq(dnaseq)
                for codon in codons:
                    df_dna[codon]=df_dna.dnaseq.apply(codon_freq, args =(codon,))
                df=pd.merge(df_prot,df_dna,on='shortid')

                # Parsing uniprot - today only Pfam
                #print "Parse uniprot"
                (dic_pfam,dic_numdoms,dic_kingdom) = parse_uniprot(uniprot_data_file)
                df["Pfam"] = df['query_id'].map(dic_pfam)
                df["NumDoms"] = df['query_id'].map(dic_numdoms)
                df["PfamType"] = df['query_id'].map(dic_kingdom)

                # export
                columns = ["query_id",  "length", "top-idp", "iupred_long", "iupred_short","iupred04_long", "iupred04_short","seg","ss_alpha", "ss_beta", "ss_coil", "ss_turn","hessa","Pfam","NumDoms","PfamType","scampi_i","scampi_o","scampi_m"]
                columns += ["freq_" + aa for aa in aas]
                columns += ["GC1","GC2","GC3","GCcoding"]
                columns += nucleotides
                columns += nucleotidepos
                columns += codons
                #for type in ["I","O","M"]:
                #    for aa in  ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
                #        colums+=type+"-"+aa

                df[columns].to_csv(out_annotation_file, index = False)

        
# load the list of PFAM domains that are shared by at least 5 bacteria and 5 euks
out_domain_ids_filename = dir+"bin/pfam_ids_orthologs_10.list"
shared_domains_pfam_ids = set(filter(None, open(out_domain_ids_filename).read().split("\n")))
shared_domains={}
for key in shared_domains_pfam_ids:
    shared_domains[key]=1



#annotate_genome('/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/data/proteomes/UP000001554_7739.fasta')
#annotate_genome('/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/data/proteomes/UP000002311_559292.fasta')
#sys.exit()
    
file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        if f.find("DNA") == -1:
            file_list += [input_dir + f]

for f in file_list:
    #print (f)
    try:
        annotate_genome(f)
    except:
        print ("ERROR on " + f)





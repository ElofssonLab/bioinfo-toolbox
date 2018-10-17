
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



dic_aa = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY',
          'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M':
          'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO',
          'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V':
          'VAL', 'Y': 'TYR'}


import string

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
       'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def get_phylum(s):
    try:
        return s.split(";")[1]
    except:
        return None    


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
        for  db in record.cross_references:
            if ( db[0] == "Pfam"):
                dic_pfam[db[1]]=1
    return dic_pfam



def extract_pfam(f,arch,bac,euk,alla):

    if not os.path.exists("./stop"):
        # check the existence of these three data files (produced separately)
        uniprot_data_file = re.sub(r'\.fasta','.txt',f)
        tempname=re.match(r'.*UP.*\_(\d+)\.fasta.*',f)
        #tax_id = int(filename.split("_")[1].split(".")[0])
        tax_id = int(tempname.group(1))

        name = taxid2name.get(tax_id,pd.np.nan)
        phylum = taxid2phylum.get(tax_id,pd.np.nan)
        kingdom = taxid2kingdom.get(tax_id,pd.np.nan)
        GC = taxid2gc.get(tax_id,pd.np.nan)
        
        (dic_pfam) = parse_uniprot(uniprot_data_file)
        f_name = f.split("/")[-1]
        for pfam in dic_pfam:
            #print ("Key:",pfam)
            alla[pfam]=1
            if (not pfam in arch.keys()):
                arch[pfam]=0
            if (not pfam in euk.keys()):
                euk[pfam]=0
            if (not pfam in bac.keys()):
                bac[pfam]=0
        if (kingdom=='Archaea'):
            for pfam in dic_pfam:
                arch[pfam]+=1
        elif (kingdom=='Bacteria'):
            for pfam in dic_pfam:
                bac[pfam]+=1
        elif (kingdom=='Eukaryota'):
            for pfam in dic_pfam:
                euk[pfam]+=1

    return(arch,bac,euk)
            
df_reference_file = data_dir  + "df_reference.csv"

if not os.path.exists(df_reference_file):

    ncbi_euks_file = data_dir + "eukaryotes.txt"
    ncbi_proks_file = data_dir + "prokaryotes.txt"
    ncbi_virs_file = data_dir + "viruses.txt"
    df_proks = pd.read_csv(ncbi_euks_file, sep = "\t")
    df_euks = pd.read_csv(ncbi_proks_file, sep = "\t")
    df_virs = pd.read_csv(ncbi_virs_file, sep = "\t")
    df_reference = pd.concat([df_euks, df_proks, df_virs])

    df_reference.to_csv(df_reference_file, index = False)

else:
    df_reference = pd.read_csv(df_reference_file)

genomes_overview_file = data_dir + "overview.txt"
df_overview = pd.read_csv(genomes_overview_file, sep="\t")
taxid2gc = df_reference.set_index("TaxID").to_dict()["GC%"]
taxid2group = df_reference[["TaxID","Group"]].set_index("TaxID").to_dict()["Group"]
group2kingdom = df_overview[["Group", "Kingdom"]].set_index("Group").to_dict()["Kingdom"]
taxid2kingdom = {}

for t in taxid2group:
    gr = taxid2group[t]
    taxid2kingdom[t] = group2kingdom[gr]


df_taxonomy = pd.read_csv(data_dir + "taxonomy.tab", sep="\t")


df_taxonomy["Phylum"] = df_taxonomy.Lineage.apply(get_phylum)
taxid2phylum = df_taxonomy.set_index("Taxon").to_dict()["Phylum"]
taxid2name = df_reference[["#Organism/Name","TaxID"]].set_index("TaxID").to_dict()['#Organism/Name']
name2taxid = df_reference[["#Organism/Name","TaxID"]].set_index("#Organism/Name").to_dict()['TaxID']

file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        if f.find("DNA") == -1:
            file_list += [input_dir + f]


arch={}
euk={}
bac={}
alla={}
#file_list=['data/proteomes/UP000001554_7739.fasta']
for f in file_list:
    try:
        extract_pfam(f,arch,bac,euk,alla)
    except:
        print ("ERROR on " + f)
cutoff=10
for pfam in sorted(alla):
    numkingdoms=0
    if (bac[pfam]>=cutoff):
        numkingdoms+=1
    if (arch[pfam]>=cutoff):
        numkingdoms+=1
    if (euk[pfam]>=cutoff):
        numkingdoms+=1
    print pfam,bac[pfam],arch[pfam],euk[pfam],numkingdoms



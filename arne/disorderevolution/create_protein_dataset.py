#/usr/bin/env python
import pandas as pd
import os
import re
import math
import numpy as np



def parse_annotation(filename):
    tempname=re.match(r'.*UP.*\_(\d+)\.fasta.*',filename)
    #tax_id = int(filename.split("_")[1].split(".")[0])
    tax_id = int(tempname.group(1))
    #print (filename,tax_id)

    df = pd.read_csv(filename)
    newdf= pd.DataFrame()
    if (len(df)==0):
        return newdf
    #df.rename(columns={'longid':'query_id'}, inplace=True)

   
    #newdf=df.filter([columns],axis=1)  # Cannot get this to work
    newdf=df[columns]  # THis gives a warning
    
    # Columns, the same for each entry
    newdf.loc[:,"taxon_id"]=tax_id
    #newdf.loc[:,'taxon_id'] = pd.Series(tax_id, index=newdf.index)
    #newdf = newdf.assign(taxon_id=pd.Series(tax_id).values)
    #newdf.loc[:,"name"] = taxid2name.get(tax_id,pd.np.nan)
    #newdf.loc[:,"phylum"] = taxid2phylum.get(tax_id,pd.np.nan)
    newdf.loc[:,"kingdom"] = taxid2kingdom.get(tax_id,pd.np.nan)
    try:
        GC = float(taxid2gc.get(tax_id,pd.np.nan))
    except:
        GC=float('nan')
    newdf.loc[:,"GC"]=GC
    try:
        GenomeSize=float(taxid2Mb.get(tax_id,pd.np.nan))*1000000
    except:
        GenomeSize=float('nan')
    newdf.loc[:,"GenomeSize"]=GenomeSize
 
    if (not (math.isnan(GC) or math.isnan(GenomeSize))):
        sumlen=df["length"]*df["GCcoding"]
        
        AveGCcoding=sumlen.sum()/df["length"].sum()
    else:
        AveGCcoding=float('nan')
    newdf.loc[:,"AveGCcoding"]=AveGCcoding
    #print (newdf)
    return newdf


def get_phylum(s):
    try:
        return s.split(";")[1]
    except:
        return None    



aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
       'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


scales =  ["ss_alpha", "ss_beta", "ss_coil", "ss_turn", "top-idp", "hessa"]

dir='/fast/arnee/proteomes/'
# Columns to read
#columns = ["length", "top-idp", "iupred_long", "iupred_short","iupred04_long", "iupred04_short","seg","ss_alpha", "ss_beta", "ss_coil", "ss_turn","hessa","scampi_i","scampi_m","scampi_o"]
columns =  ["query_id","length"]
columns += ["freq_" + aa for aa in aas]
columns += ["GC1","GC2","GC3","GCcoding"]

if  (not os.path.isdir(dir)):
    dir='/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/'

data_dir = dir + "/data/"
results_dir = dir + "/results/"
input_dir = results_dir + "extended/"



# create and read the reference table file (to extract GC% and other metrics for each species)
df_reference_file = data_dir  + "df_reference.csv"

if not os.path.exists(df_reference_file):

    ncbi_euks_file = data_dir + "eukaryotes.txt"
    ncbi_proks_file = data_dir + "prokaryotes.txt"
    ncbi_virs_file = data_dir + "viruses.txt"
    genomes_overview_file = data_dir + "overview.txt"

    if not os.path.exists(ncbi_euks_file):
        ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"
        cmd = "wget '" + ncbi_url + "' -O " + ncbi_euks_file
        os.system(cmd)

    if not os.path.exists(ncbi_proks_file):
        ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
        cmd = "wget '" + ncbi_url + "' -O " + ncbi_proks_file
        os.system(cmd)

    if not os.path.exists(ncbi_virs_file):
        ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/viruses.txt"
        cmd = "wget '" + ncbi_url + "' -O " + ncbi_virs_file
        os.system(cmd)

    df_proks = pd.read_csv(ncbi_euks_file, sep = "\t")
    df_euks = pd.read_csv(ncbi_proks_file, sep = "\t")
    df_virs = pd.read_csv(ncbi_virs_file, sep = "\t")
    df_reference = pd.concat([df_euks, df_proks, df_virs])

    
else:
    df_reference = pd.read_csv(df_reference_file)

genomes_overview_file = data_dir + "overview.txt"
if not os.path.exists(genomes_overview_file):
    ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"
    cmd = "wget '" + ncbi_url + "' -O " + genomes_overview_file

df_overview = pd.read_csv(genomes_overview_file, sep="\t")



taxid2gc = df_reference.set_index("TaxID").to_dict()["GC%"]
taxid2Kb = df_reference.set_index("TaxID").to_dict()["Size (Kb)"]
taxid2Mb = df_reference.set_index("TaxID").to_dict()["Size (Mb)"]



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
    if f.endswith("_annotation.csv"):
        file_list += [input_dir + f]
            
#file_list = [input_dir+"UP000005777_641146.fasta_annotation.csv",input_dir+"UP000000466_1117647.fasta_annotation.csv"]


#allcolumns=columns+["taxon_id","name","phylum","kingdom","GC","GenomeSize","AveGCcoding"]
data=pd.DataFrame()

for f in file_list:
    print (f)
    d = parse_annotation(f)
    if len(data)>0:
        #print (d)
        #print (data)
        data=pd.concat([data,d])
        #print (data)
    else:
        data=d
    #print (data)
data.to_csv(results_dir + "df_GC-freq.csv", index = False)


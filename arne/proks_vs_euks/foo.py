#!/usr/bin/env python
import pandas as pd
import os



aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

data_dir = "../data/"
results_dir = "../results/"
input_dir = results_dir + "/uniprot_pfam_single/"
#input_dir = results_dir + "/test/"

file_list = []
for f in os.listdir(input_dir):
    if f.endswith("_full.csv"):
        file_list += [input_dir + f]
            
print file_list



columns = ["query_id",  "length_translation", "top-idp", "hessa", "iupred_long", "iupred_short", "seg"]
columns += ["freq_" + aa for aa in aas]
columns += ["ss_alpha", "ss_beta", "ss_coil", "ss_turn"]
columns += ["pfamA_acc"]

# Add the new columns
newcolumns = list(columns)
newcolumns += ["tax_id","phylum","kingdom","GC"]

data = []
for f in file_list:
    #d = parse_annotation(f)
    print "F",f
    df = pd.read_csv(f,names=columns)
    print "DF",df
    n_proteins=len(df.query_id)
    for index, row in df.iterrows():
        print "ROW",row
        mtemp=row["query_id"]
        print "MTEMP",mtemp
        mnemonic=mtemp.split("_")[1]
        print "MNEM",mnemonic
        tax_id=mnemonic
        phylum = mnemonic
        kingdom = mnemonic
        GC = 25.
        row["tax_id"]=tax_id
        row["phylum"]=phylum
        row["kingdom"]=kingdom
        row["GC"]=GC
        #row+=[tax_id,phylum,kingdom,GC]
        data += [row]


print data
df_final = pd.DataFrame(data,columns=newcolumns)
df_final.to_csv("foo.csv", index = False,header = True )



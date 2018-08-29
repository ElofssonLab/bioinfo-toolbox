#!/usr/bin/env python
import pandas as pd
import os


def get_phylum(s):
    try:
        return s.split(";")[1]
    except:
        return None    

def get_kingdom(s):
    try:
        return s.split(";")[0]
    except:
        return None    

def parse_annotation(filename):

    tax_id = int(filename.split("_")[1].split(".")[0])
    print tax_id

    df = pd.read_csv(filename)
    n_proteins = len(df.query_id)

    df_mean = df.mean()
    
    columns = ["length", "top-idp", "iupred_long", "iupred_short","seg","ss_alpha", "ss_beta", "ss_coil", "ss_turn","hessa"]
    columns += ["freq_" + aa for aa in aas]
    columns += ["pfamA_acc"]
    
    ret_dic = {}    
    ret_dic["taxon_id"] = tax_id
    ret_dic["phylum"] = taxid2phylum.get(tax_id,pd.np.nan)
    ret_dic["kingdom"] = taxid2kingdom.get(tax_id,pd.np.nan)
    ret_dic["GC"] = taxid2gc.get(tax_id,pd.np.nan)
    ret_dic["count_protein"] = n_proteins

    for c in columns:
        ret_dic[c] = df_mean[c]
    
    #gcs = df_reference.loc[df_reference["TaxID"] == tax_id]["GC%"].astype(float)
    #ret_dic["GC"] = pd.np.mean(list(gcs))

    return ret_dic


############################################################################
##    MAIN
############################################################################

    

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

data_dir = "../data/"
results_dir = "../results/"
input_dir = results_dir + "/uniprot_pfam_single/"
#input_dir = results_dir + "/test/"


# create and read the reference table file (to extract GC% and other metrics for each species)
df_reference_file = data_dir  + "df_reference.csv"

if not os.path.exists(df_reference_file):
    print "Getting reference File"

    ncbi_euks_file = data_dir + "eukaryotes.txt"
    ncbi_proks_file = data_dir + "prokaryotes.txt"
    ncbi_virs_file = data_dir + "viruses.txt"
   
    if not os.path.exists(ncbi_euks_file):
        print "Getting eukaryotic file"
        ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"
        cmd = "wget '" + ncbi_url + "' -O " + ncbi_euks_file
        os.system(cmd)

    if not os.path.exists(ncbi_proks_file):
        print "Getting prokaryotic file"
        ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
        cmd = "wget '" + ncbi_url + "' -O " + ncbi_proks_file
        os.system(cmd)

    if not os.path.exists(ncbi_virs_file):
        print "Getting Virus file"
        ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/viruses.txt"
        cmd = "wget '" + ncbi_url + "' -O " + ncbi_virs_file
        os.system(cmd)

    df_proks = pd.read_csv(ncbi_euks_file, sep = "\t")
    df_euks = pd.read_csv(ncbi_proks_file, sep = "\t")
    df_virs = pd.read_csv(ncbi_virs_file, sep = "\t")
    df_reference = pd.concat([df_euks, df_proks, df_virs])

    df_reference.to_csv(df_reference_file, index = False)

else:
    print "Reading Reference File",df_reference_file
    df_reference = pd.read_csv(df_reference_file)

genomes_overview_file = data_dir + "overview.txt"
if not os.path.exists(genomes_overview_file):
    print "getting overiew file"
    ncbi_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"
    cmd = "wget '" + ncbi_url + "' -O " + genomes_overview_file


print "reading overview file",genomes_overview_file
df_overview = pd.read_csv(genomes_overview_file, sep="\t")
print "done"


taxid2gc = df_reference.set_index("TaxID").to_dict()["GC%"]
#print taxid2gc
#print "1117:",taxid2gc["1117"]
#print "27334:",taxid2gc["sl334"]
taxid2group = df_reference[["TaxID","Group"]].set_index("TaxID").to_dict()["Group"]
group2kingdom = df_overview[["Group", "Kingdom"]].set_index("Group").to_dict()["Kingdom"]
taxid2kingdom = {}

#print "taxome classification"
#for t in taxid2group:
#    gr = taxid2group[t]
#    taxid2kingdom[t] = group2kingdom[gr]


print "reading taxonomy tab"
df_taxonomy = pd.read_csv(data_dir + "taxonomy.tab", sep="\t")
print "Set taxon index"
memo2taxid = df_taxonomy.set_index("Mnemonic").to_dict()["Taxon"]
print "Get phylym"
df_taxonomy["Phylum"] = df_taxonomy.Lineage.apply(get_phylum)
print "Setting phylum index"
taxid2phylum = df_taxonomy.set_index("Taxon").to_dict()["Phylum"]

df_taxonomy["Kingdom"] = df_taxonomy.Lineage.apply(get_kingdom)
print "Setting phylum index"
taxid2kingdom = df_taxonomy.set_index("Taxon").to_dict()["Kingdom"]


datasets = ["full",  "domains_shared", "linkers"]
columns = ["query_id",  "length_translation", "top-idp", "hessa", "iupred_long", "iupred_short", "seg"]
columns += ["freq_" + aa for aa in aas]
columns += ["ss_alpha", "ss_beta", "ss_coil", "ss_turn"]
columns += ["pfamA_acc"]

# Add the new columns
newcolumns = list(columns)
newcolumns += ["tax_id","phylum","kingdom","GC"]

for set in datasets:
    print "getting input files" 
    file_list = []
    for f in os.listdir(input_dir):
        if f.endswith("_"+set+".csv"):
            file_list += [input_dir + f]
    data = []
    for f in file_list:
        print "reading ",f
        #d = parse_annotation(f)
        df = pd.read_csv(f,names=columns)
        n_proteins=len(df.query_id)
    
        for index, row in df.iterrows():
            mtemp=row["query_id"]
            mnemonic=mtemp.split("_")[1]
            tax_id=memo2taxid.get(mnemonic,pd.np.nan)
            phylum = taxid2phylum.get(tax_id,pd.np.nan)
            kingdom = taxid2kingdom.get(tax_id,pd.np.nan)
            GC = taxid2gc.get(tax_id,pd.np.nan)
            row["tax_id"]=tax_id
            row["phylum"]=phylum
            row["kingdom"]=kingdom
            row["GC"]=GC
            #row+=[tax_id,phylum,kingdom,GC]
            data += [row]


    df_final = pd.DataFrame(data,columns=newcolumns)
    df_final.to_csv(results_dir + "df_pfam_single_proteomes_with_pfam_and_kingdom_"+set+".csv", index = False,header = True )



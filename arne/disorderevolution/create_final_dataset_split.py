import pandas as pd
import os
import re


def parse_annotation(filename,ty):
    tempname=re.match(r'.*UP.*\_(\d+)\.fasta.*',filename)
    #tax_id = int(filename.split("_")[1].split(".")[0])
    tax_id = int(tempname.group(1))
    print (filename,tax_id,ty)

    fulldf = pd.read_csv(filename)
    tempdf=fulldf[fulldf['PfamType'].notnull()]
    if (len(tempdf)==0):
        ret_dic = {}    
        sum_dic = {}    
        return ret_dic,sum_dic
    if (ty == "All"):
        df = tempdf.copy()
    else:
        df = tempdf.loc[(tempdf.PfamType == ty) ]
    if (len(df)==0):
        ret_dic = {}    
        sum_dic = {}    
        return ret_dic,sum_dic
    n_proteins = len(df.query_id)
    
    df_sum = df.sum()
    df_mean = df.mean()
    
    columns = ["length", "top-idp", "iupred_long", "iupred_short","seg","ss_alpha", "ss_beta", "ss_coil", "ss_turn","hessa"]
    columns += ["freq_" + aa for aa in aas]

    ret_dic = {}    
    ret_dic["taxon_id"] = tax_id
    ret_dic["name"] = taxid2name.get(tax_id,pd.np.nan)
    ret_dic["phylum"] = taxid2phylum.get(tax_id,pd.np.nan)
    ret_dic["kingdom"] = taxid2kingdom.get(tax_id,pd.np.nan)
    ret_dic["GC"] = taxid2gc.get(tax_id,pd.np.nan)
    ret_dic["GC_genomic"] = taxid2gc.get(tax_id,pd.np.nan)
    ret_dic["count_protein"] = n_proteins
    ret_dic["length_translation"]=df["length"].sum()
    sum_dic = {}    
    sum_dic["taxon_id"] = tax_id
    sum_dic["name"] = taxid2name.get(tax_id,pd.np.nan)
    sum_dic["phylum"] = taxid2phylum.get(tax_id,pd.np.nan)
    sum_dic["kingdom"] = taxid2kingdom.get(tax_id,pd.np.nan)
    sum_dic["GC"] = taxid2gc.get(tax_id,pd.np.nan)
    sum_dic["GC_genomic"] = taxid2gc.get(tax_id,pd.np.nan)
    sum_dic["count_protein"] = n_proteins
    sum_dic["length_translation"]=df["length"].sum()

    for c in columns:
        #df.loc[:,c+"-sum"] = df.loc[:,c]*df.loc[:,"length"]*n_proteins
        if (c=="length"):
            df[c+"-sum"] = df[c]
        elif ( c == "length_translation"):
            df[c+"-sum"] = df[c]            
        #elif ( c in scales ):
        #    df[c+"-sum"] = df[c]
        else:
            df[c+"-sum"] = df[c]*df["length"]

        if (c=="length"):
            sum_dic[c] = df[c+"-sum"].mean()
        else:
            sum_dic[c] = df[c+"-sum"].sum()
        ret_dic[c] = df_mean[c]
        
    #gcs = df_reference.loc[df_reference["TaxID"] == tax_id]["GC%"].astype(float)
    #ret_dic["GC"] = pd.np.mean(list(gcs))

    return ret_dic,sum_dic


def get_phylum(s):
    try:
        return s.split(";")[1]
    except:
        return None    



aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


scales =  ["ss_alpha", "ss_beta", "ss_coil", "ss_turn", "top-idp", "hessa"]

dir='/scratch2/arne/annotate_uniprot_proteomes/'

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
            
#file_list = [input_dir+"UP000077428_66851.fasta_annotation.csv"]



for ty in ["All","Shared","None","Unique"]:
    data = []
    summ = []
    for f in file_list:
        d,s = parse_annotation(f,ty)
        data += [d]
        summ += [s]
    df_final = pd.DataFrame(data)
    df_final.to_csv(results_dir + "df_uniprot_reference_proteomes_per_species-"+ty+".csv", index = False)
    df_sum = pd.DataFrame(summ)
    df_sum.to_csv(results_dir + "df_uniprot_reference_proteomes_per_species-sum-"+ty+".csv", index = False)


import pandas as pd
import os
import re

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']



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

    df_reference.to_csv(df_reference_file, index = False)

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

def get_phylum(s):
    try:
        return s.split(";")[1]
    except:
        return None    

df_taxonomy["Phylum"] = df_taxonomy.Lineage.apply(get_phylum)
taxid2phylum = df_taxonomy.set_index("Taxon").to_dict()["Phylum"]
taxid2name = df_reference[["#Organism/Name","TaxID"]].set_index("TaxID").to_dict()['#Organism/Name']
name2taxid = df_reference[["#Organism/Name","TaxID"]].set_index("#Organism/Name").to_dict()['TaxID']




file_list = []
for f in os.listdir(input_dir):
    if f.endswith("_annotation.csv"):
        file_list += [input_dir + f]
            


def parse_annotation(filename):
    tempname=re.match(r'.*UP.*\_(\d+)\.fasta.*',filename)
    #tax_id = int(filename.split("_")[1].split(".")[0])
    tax_id = int(tempname.group(1))
    print (filename,tax_id)

    df = pd.read_csv(filename)
    n_proteins = len(df.query_id)

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

    for c in columns:
        ret_dic[c] = df[c]
        ret_dic[c+"-avg"] = df_mean[c]
    
    #gcs = df_reference.loc[df_reference["TaxID"] == tax_id]["GC%"].astype(float)
    #ret_dic["GC"] = pd.np.mean(list(gcs))

    return ret_dic


data = []
for f in file_list:
    d = parse_annotation(f)
    data += [d]

df_final = pd.DataFrame(data)
df_final.to_csv(results_dir + "df_uniprot_reference_proteomes_per_species-All.csv", index = False)



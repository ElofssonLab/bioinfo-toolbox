import os
import pandas as pd

# create the relevant directories
data_dir = "../data/"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

proteomes_dir = data_dir + "proteomes/"
if not os.path.exists(proteomes_dir):
    os.makedirs(proteomes_dir)


# download the list of reference proteomes from uniprot
outfile_list = data_dir + 'uniprot_reference_proteomes.tsv'

if not os.path.exists(outfile_list):
    list_url = 'http://www.uniprot.org/proteomes/?sort=&desc=&compress=no&query=&fil=reference:yes&force=no&preview=false&format=tab&columns=id,name,organism-id,proteincount'

    cmd = "wget '" + list_url + "' -O " + outfile_list
    os.system(cmd)




ref_proteomes_file_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2017_11.tar.gz'
ref_proteomes_file_compressed = data_dir + "Reference_Proteomes_2017_11.tar.gz"


if not os.path.exists(ref_proteomes_file_compressed):
    cmd = "wget '" + ref_proteomes_file_url+ "' -O " + ref_proteomes_file_compressed
    os.system(cmd)



'''
# read the csv file with the list
df = pd.read_csv(outfile_list, sep="\t")

for i,r in df.iterrows():
    taxid = str(r["Organism ID"])
    outfname = proteomes_dir + taxid + "_proteome.fasta"
    if not os.path.exists(outfname):
        p_name = str(r["Proteome ID"])

        proteome_url = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=proteome:'+p_name+'&fil=&force=no&preview=false&format=fasta'

        cmd = "wget '" + proteome_url + "' -O " + outfname
        os.system(cmd)
'''

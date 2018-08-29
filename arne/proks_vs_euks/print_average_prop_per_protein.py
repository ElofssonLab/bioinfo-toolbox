import pandas as pd

prop = "iupred_long"
ty = "domains_orthologs_5"
df = pd.read_feather("../results/final_datasets/pfam-" + ty + "_annotations.feather")[[prop,"kingdom"]].dropna()


## eliminate bad taxa
genomes_overview_file = data_dir + "overview.txt"
genomes_euks_file = data_dir + "eukaryotes.txt"
genomes_proks_file = data_dir + "prokaryotes.txt"

df_overview = pd.read_csv(genomes_overview_file, sep="\t")
df_proks = pd.read_csv(genomes_proks_file, sep = "\t")
df_euks = pd.read_csv(genomes_euks_file, sep = "\t")

df_reference = pd.concat([df_euks, df_proks])

name2taxid = df_reference[["#Organism/Name","TaxID"]].set_index("#Organism/Name").to_dict()['TaxID']

# these taxa should be excluded because their genetic code is different, so translation is wrong
mycotaxa = []
for t in name2taxid:
	tax = t.split()[0]
	if tax in ["Mycoplasma","Spiroplasma","Ureaplasma","Mesoplasma"]:
		mycotaxa += [name2taxid[t]]

good_taxa = set(df.taxon_id.unique()) - set(mycotaxa)

df = df.loc[df.taxon_id.isin(good_taxa)]
##


print ty
for k in ["Eukaryotes", "Bacteria"]:
	print k, df.loc[df.kingdom == k][prop].mean()
	


#print ty
#for k in ["Eukaryotes", "Bacteria"]:
#	print k, df.loc[df.kingdom == k][prop].mean()
	

import pandas as pd

good_taxa = filter(None, open("../data/pfam/good_taxa").read().split("\n"))

prop = "iupred_long"

datasets = ["full", "alignments", "domains", "linkers"]
datasets += ["domains_exclusive_euk","domains_exclusive_bac","domains_orthologs_5"]


for ty in datasets:

    infname = "../results/final_datasets/pfam-" + ty + "_annotations.feather"
    outfname = "../results/final_datasets/derived/" + prop + "_" + ty + ".csv"
    print prop
    print "IN:", infname
    print "OUT:", outfname
    print

    df = pd.read_feather(infname)[[prop,"taxon_id","kingdom"]].dropna()

    df.taxon_id = df.taxon_id.astype(int)

    df = df.loc[df.taxon_id.isin(good_taxa)]

    df.reset_index().to_csv(outfname, index=False)






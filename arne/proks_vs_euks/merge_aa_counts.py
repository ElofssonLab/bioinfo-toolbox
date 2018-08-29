
import pandas as pd
import os


ty = "linkers"


if not os.path.exists("../results/final_datasets/pfam-"+ty+"_annotations.csv"):
    df = pd.read_feather("../results/final_datasets/tmp/pfam-"+ty+"_annotations.feather")
    print "opened main df"

    df_counts = pd.read_feather("../results/final_datasets/tmp/pfam-"+ty+"_annotations_aa_counts.feather")
    print "opened count df"

    df = pd.merge(df, df_counts, on="query_id")
    print "merged"

    df.to_csv("../results/final_datasets/pfam-"+ty+"_annotations.csv")
    print "exported to CSV"

else:

    df = pd.read_csv("../results/final_datasets/pfam-"+ty+"_annotations.csv")

    df.reset_index().to_feather("../results/final_datasets/pfam-"+ty+"_annotations.feather")
    print "exported to feather"



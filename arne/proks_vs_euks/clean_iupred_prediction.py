
import os
import commands



input_dir = "../results/uniprot_sequences_groups/"
annotations_dir = "../results/uniprot_pfam_annotations/"


for ty in ["long", "short"]:
    file_list = []
    for f in os.listdir(input_dir):
        if f.endswith("_iupred_" + ty):
            file_list += [input_dir + f]
       
    for f in file_list:
        # clean empty files
        st = os.stat(f)
        if st.st_size == 0:
            os.system("rm " + f)

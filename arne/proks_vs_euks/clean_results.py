
import os
import commands



input_dir = "../results/uniprot_sequences_groups/"
annotations_dir = "../results/uniprot_pfam_annotations/"

file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        file_list += [input_dir + f]
        
for f in file_list:
    
    f_name = f.split("/")[-1]
    out_annotation_file = annotations_dir + f_name +"_annotation_"

    #print f, out_annotation_file
    
    if os.path.exists(out_annotation_file):
        csv_fname_full = out_annotation_file + "full.csv"
        csv_fname_linkers = out_annotation_file + "linkers.csv"
        csv_fname_domains_shared = out_annotation_file + "domains_shared.csv"
        csv_fname_domains_other = out_annotation_file + "domains_other.csv"
        
        res = os.path.isfile(csv_fname_full) & os.path.isfile(csv_fname_linkers) & os.path.isfile(csv_fname_domains_shared) & os.path.isfile(csv_fname_domains_other)
        
        if res == False:
            print f
            for outf in [out_annotation_file, csv_fname_full,csv_fname_linkers, csv_fname_domains_other, csv_fname_domains_shared]:
                cmd = "rm " +outf
                #os.system(cmd)

#!/usr/bin/env python

import os
import commands



input_dir = "/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes//data/proteomes/"

scampi= "/pfs/nobackup/home/a/arnee/git/scampi2/bin/scampi/SCAMPI_run.pl"

#~ file_list = []
#~ for f in os.listdir(input_dir):
    #~ if f.endswith("_iupred_short"):
        #~ file_list += [input_dir + f]
   
#~ for f in file_list:
    #~ cmd = "mv " + f + " " + f.replace("_iupred_short", "_iupred_short_raw")
    #~ print cmd
    #~ os.system(cmd)
   
file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        if f.find("DNA") == -1:
            file_list += [input_dir + f]
        

for f in file_list:
    
    scampi_data_file = f +  ".scampi" 
        
    # clean empty files
    if os.path.exists(scampi_data_file):
        st = os.stat(scampi_data_file)
        #if st.st_size == 0:
        #    os.system("rm " + scampi_data_file)

    if not os.path.exists(scampi_data_file):
        print "creating " + scampi_data_file
        os.system("touch " + scampi_data_file)
        os.system(scampi + "  " + f + " " + scampi_data_file)

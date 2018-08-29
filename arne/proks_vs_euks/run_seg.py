
import os
import commands


############################################################################
##    MAIN
############################################################################


input_dir = "../results/uniprot_sequences_groups/"
#annotations_dir = "../results/uniprot_pfam_annotations/"

    
file_list = []
for f in os.listdir(input_dir):
    if f.endswith(".fasta"):
        file_list += [input_dir + f]
        

for f in file_list:
    print "processing " + f
                       
    # SEG, for the whole, full sequences file
    seg_data_raw_file = f +  ".data_seg"

    if os.path.exists(seg_data_raw_file):
        st = os.stat(seg_data_raw_file)
        if st.st_size == 0:
            os.system("rm " + seg_data_raw_file)

    if not os.path.exists(seg_data_raw_file):
        cmd = './seg '+ f + " -x > " + seg_data_raw_file
        
        print cmd
        try:
            ret = commands.getoutput(cmd)
            if ret != "":
                print "ERROR: ", ret

        except:
            print "ERROR: could not run SEG (probably could not allocate memory)"


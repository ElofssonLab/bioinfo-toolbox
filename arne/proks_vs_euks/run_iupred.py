
import os
import commands



input_dir = "/pfs/nobackup/home/w/wbasile/proks_euks/results/uniprot_sequences_groups/"

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
        file_list += [input_dir + f]
        

for f in file_list:
    #print f
    for iupred_param in ["long", "short"]:
        iupred_data_file = f +  ".data_iupred_" + iupred_param 
        iupred_data_file_04 = f +  ".data_iupred_0.4_" + iupred_param 
        iupred_data_file_raw = f +  ".data_iupred_"+iupred_param+"_raw" 
        
        # This makes it possible to run in parallell
        if os.path.exists(iupred_data_file) and  os.path.exists(iupred_data_file_04):
            #print "Skipping ", f
            continue
        # clean empty files
        #if os.path.exists(iupred_data_file):
        #    st = os.stat(iupred_data_file)
        #    if st.st_size == 0:
        #        os.system("rm " + iupred_data_file)
        #
        #if os.path.exists(iupred_data_file_04):
        #    st = os.stat(iupred_data_file_04)
        #    if st.st_size == 0:
        #        os.system("rm " + iupred_data_file_04)


        if not os.path.exists(iupred_data_file):
            
            print "creating " + iupred_data_file
            
            os.system("touch " + iupred_data_file)
        
            if os.path.exists(iupred_data_file_raw):

                # remove empty files
                st = os.stat(iupred_data_file_raw)
                if st.st_size == 0:
                    print "Removing file, size 0",iupred_data_file_raw
                    os.system("rm " + iupred_data_file_raw)
                else:
                    group_num = int(f.split(".fasta")[0].split("_")[-1])
                    
                    # remove the incomplete raw files
                    if group_num != 1208:
                        cmd = "grep -c Prediction " + iupred_data_file_raw 
                        print cmd
                        ret = int(commands.getoutput(cmd))
                        print "Testing number of entries:", ret 
                        if ret != 50000:
                            print ("rm " + iupred_data_file_raw)
                            #sys.exit()
                            os.system("rm " + iupred_data_file_raw)

            if not os.path.exists(iupred_data_file_raw):
            
                cmd = "./iupred_multi " + f + " " + iupred_param + " > " + iupred_data_file_raw
                print cmd
                os.system(cmd)

                
            
                
            #PARSE
            print "parsing " + iupred_data_file_raw + " into " + iupred_data_file
            ps = filter(None, open(iupred_data_file_raw).read().split("# Prediction output "))

            iupred_dic = {}
            for p in ps:
                #print p
                query_id = ""
                values = []

                lines = filter(None,p.split("\n"))
                for line in lines:
                    if line[0] == "#":
                        query_id = line.split()[1]
                    else:
                        values += [float(line.split()[-1])]

                # calculate the disorder as the fraction of disordered (>0.5) residues
                if query_id != "":
                    
                    diso_aa = ""
                    for a in values:
                        if a > 0.5:
                            diso_aa += "x"
                        else:
                            diso_aa += "n"
                    
                    iupred_dic[query_id] = diso_aa
                    
                    
            with open(iupred_data_file, "w") as outf:
                for k in iupred_dic.keys():
                    outf.write(">" + k + "\n" + iupred_dic[k] + "\n")

        if not os.path.exists(iupred_data_file_04):
            print "creating " + iupred_data_file_04
            os.system("touch " + iupred_data_file_04)
            if os.path.exists(iupred_data_file_raw):

                # remove empty files
                st = os.stat(iupred_data_file_raw)
                if st.st_size == 0:
                    print "Removing file, size 0",iupred_data_file_raw
                    os.system("rm " + iupred_data_file_raw)
                else:
                    group_num = int(f.split(".fasta")[0].split("_")[-1])
                    
                    # remove the incomplete raw files
                    if group_num != 1208:
                        cmd = "grep -c Prediction " + iupred_data_file_raw
                        ret = int(commands.getoutput(cmd))
                        
                        if ret != 50000:
                               print ("rm " + iupred_data_file_raw)
                               #sys.exit()
                               os.system("rm " + iupred_data_file_raw)

            if not os.path.exists(iupred_data_file_raw):
            
                cmd = "./iupred_multi " + f + " " + iupred_param + " > " + iupred_data_file_raw
                print cmd
                os.system(cmd)
        
                            
                        
            #PARSE
            print "parsing " + iupred_data_file_raw  + " into " + iupred_data_file
            ps = filter(None, open(iupred_data_file_raw).read().split("# Prediction output "))

            iupred_dic = {}
            for p in ps:
                query_id = ""
                values = []

                lines = filter(None,p.split("\n"))
                for line in lines:
                    if line[0] == "#":
                        query_id = line.split()[1]
                    else:
                        values += [float(line.split()[-1])]

                # calculate the disorder as the fraction of disordered (>0.5) residues
                if query_id != "":
                    
                    diso_aa = ""
                    for a in values:
                        if a > 0.4:
                            diso_aa += "x"
                        else:
                            diso_aa += "n"
                    
                    iupred_dic[query_id] = diso_aa
                    
                    
            with open(iupred_data_file_04, "w") as outf:
                for k in iupred_dic.keys():
                    outf.write(">" + k + "\n" + iupred_dic[k] + "\n")
        #sys.exit()
            


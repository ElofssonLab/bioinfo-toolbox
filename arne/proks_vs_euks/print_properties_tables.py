
import pickle
import os

out_dic_fname = "props_dic.pickle"

props = ["iupred_long","iupred_short", "seg"]

if not os.path.exists(out_dic_fname):

    import pandas as pd

    dic_res = {}

    for ty in ["full","pfam","domains","linkers"]:

        print ty
        
        df = pd.read_feather("./results/df_" + ty + "_annotation.feather")[["kingdom","length_translation"] + props]

        dic_res[ty] = {}

        for kingdom in ["Eukaryota", "Bacteria"]:
            dic_res[ty][kingdom] = {}
            kdf = df.loc[df.kingdom == kingdom]
            len_tot = kdf["length_translation"].sum()

            dic_res[ty][kingdom]["total_residues"] = len_tot    

            for prop in props:
                print prop
                prop_tot = kdf[prop].sum()
                dic_res[ty][kingdom][prop] = prop_tot
            
    with open(out_dic_fname, 'wb') as handle:
        pickle.dump(dic_res, handle, protocol=pickle.HIGHEST_PROTOCOL)



else:

    with open(out_dic_fname, 'rb') as handle:
        dic_res = pickle.load(handle)


    for prop in props:
        print prop
        
        for kingdom in ["Eukaryota", "Bacteria"]:
            print kingdom
            #print "Total residues:", dic_res["full"][kingdom]["total_residues"] 
            print "Dataset,", ",".join([prop+" residues", "Total residues", prop+"% (local)", prop+"% (global)"])
            for ty in ["full","pfam","domains","linkers"]:
                
                perc_glob = float(dic_res[ty][kingdom][prop])/float(dic_res["full"][kingdom]["total_residues"])*100.0
                perc_loc = float(dic_res[ty][kingdom][prop])/float(dic_res[ty][kingdom]["total_residues"])*100.0

                print ty.upper() ,",", int(dic_res[ty][kingdom][prop]),",", int(dic_res[ty][kingdom]["total_residues"]),",", round(perc_loc, 2),",", round(perc_glob,2)

            print                
        print
                


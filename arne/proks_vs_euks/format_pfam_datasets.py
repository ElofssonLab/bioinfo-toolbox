# coding: utf-8
import pandas as pd
import os
import sys
import os.path

def get_kingdom(s):
    try:
        return s.split(";")[0]
    except:
        return None



def get_phylum(s):
    try:
        return s.split(";")[1]
    except:
        return None    

def set_query_id(q):
    return q.split("/")[0]


dir='/scratch2/arne/proks-vs-euk/'

if  (not os.path.isdir(dir)):
    dir='/pfs/nobackup/home/w/wbasile/proks_euks/'

data_dir = dir+"./data/"
result_dir = dir + "./results/"


# download and open the three reference files
genomes_overview_file = data_dir + "overview.txt"
genomes_euks_file = data_dir + "eukaryotes.txt"
genomes_proks_file = data_dir + "prokaryotes.txt"

df_overview = pd.read_csv(genomes_overview_file, sep="\t")

df_proks = pd.read_csv(genomes_proks_file, sep = "\t")
df_euks = pd.read_csv(genomes_euks_file, sep = "\t")
df_reference = pd.concat([df_euks, df_proks])

taxid2gc = df_reference.set_index("TaxID").to_dict()["GC%"]

df_taxonomy = pd.read_csv(data_dir+"/taxonomy.tab", sep="\t")
df_taxonomy["Kingdom"] = df_taxonomy.Lineage.apply(get_kingdom)
taxid2kingdom = df_taxonomy.set_index("Taxon").to_dict()["Kingdom"]


taxid2name = df_taxonomy.set_index("Taxon").to_dict()["Scientific name"]





df_taxonomy["Phylum"] = df_taxonomy.Lineage.apply(get_phylum)
taxid2phylum = df_taxonomy.set_index("Taxon").to_dict()["Phylum"]

for ty in  ["full", "pfam", "domains", "linkers", "N-linkers", "M-linkers", "C-linkers"]:
    #["full", "pfam", "domains", "linkers"]:

    

    ############################################################
    ######################################################3
    #ty = "full"
    #ty = "pfam"
    #ty = "domains"
    #ty = "linkers"
    ###########################################################
    ##########################################################3

    file1=result_dir+"./df_" + ty + ".csv"
    file2=result_dir+"./df_" + ty + "_annotation_with_taxid.csv"
    if os.path.getmtime(file1) < os.path.getmtime(file2):
        phase = 2
    else:
        phase = 1

    print "Phase:",phase
    ###########################################################
    ##########################################################3


    if phase == 1:

        df = pd.read_csv(result_dir+"./df_" + ty + ".csv")


        # In[4]:
        '''
        # only for the 3 new domains sets

        df["query_id_tmp"] = df.query_id.apply(set_query_id)
        df["query_domain_id"] = df.query_id
        df["query_id"] = df.query_id_tmp
        #################################
        '''

        # In[ ]:

        # annotate each sequence with its taxon id
        if "taxon_id" not in df.columns:

            df_uniprot_sel = pd.read_feather(data_dir+"/uniprot2taxid_2_fields.feather")

            #df_uniprot_sel.columns = ["index","query_id","taxon_id"]
            df = pd.merge(df,df_uniprot_sel,on="query_id",how="left")

            df = df[['query_id', 'length_translation', 'top-idp', 'hessa',
               'iupred_long', 'iupred_short', 'seg', 'freq_A', 'freq_C',
               'freq_D', 'freq_E', 'freq_F', 'freq_G', 'freq_H', 'freq_I',
               'freq_K', 'freq_L', 'freq_M', 'freq_N', 'freq_P', 'freq_Q',
               'freq_R', 'freq_S', 'freq_T', 'freq_V', 'freq_W', 'freq_Y',
               'ss_alpha', 'ss_beta', 'ss_coil',
               'ss_turn','taxon_id']]

            #df.to_csv("../data/pfam/pfam-full_annotations.csv", index=False)


        df = df.dropna(subset=["taxon_id"])

        df.taxon_id = df.taxon_id.astype(int)

        df["kingdom"] = df.taxon_id.map(taxid2kingdom)

        df["phylum"] = df.taxon_id.map(taxid2phylum)

        df["GC_genomic"] = df.taxon_id.map(taxid2gc)
        df.GC_genomic = df.GC_genomic.replace('-',pd.np.nan).astype(float)


        df.to_csv(result_dir+"./df_" + ty + "_annotation_with_taxid.csv", index=False)


    else:
        if not os.path.exists(result_dir+"./df_" + ty + "_annotation.feather"):
            df = pd.read_csv(result_dir+"./df_" + ty + "_annotation_with_taxid.csv")
            df.GC_genomic = df.GC_genomic.replace('-',pd.np.nan).astype(float)

            #df = df.dropna(subset=["iupred_long"]).reset_index()

            df.loc[df.length_translation > 0].reset_index().to_feather(result_dir+"./df_" + ty + "_annotation.feather")

        if not os.path.exists(result_dir+"./df_" + ty + "_annotation_by_species.feather"):

            if not os.path.exists(result_dir+"./df_" + ty + "_annotation_by_species.csv"):

                df = pd.read_feather(result_dir+"./df_" + ty + "_annotation.feather")

                for p in ['top-idp', 'hessa','ss_alpha', 'ss_beta', 'ss_coil',
                       'ss_turn']:

                    # multiply back these averages by the length, to obtain the raw n. of residues
                    df[p] = df[p] * df["length_translation"]


                gdf = pd.DataFrame(df.groupby("taxon_id").sum())
                gdf["count_protein"]=df.groupby("taxon_id")["taxon_id"].transform('count')
                gdf["taxon_id"] = list(gdf.index)
                gdf["kingdom"] = gdf.taxon_id.map(taxid2kingdom)
                gdf["phylum"] = gdf.taxon_id.map(taxid2phylum)
                gdf["name"] = gdf.taxon_id.map(taxid2name)

                gdf["GC_genomic"] = gdf.taxon_id.map(taxid2gc)
                gdf.GC_genomic = gdf.GC_genomic.replace('-',pd.np.nan).astype(float)

                cols = ['length_translation', 'top-idp', 'hessa',
                       'iupred_long', 'iupred_short', 'seg', 'freq_A', 'freq_C',
                       'freq_D', 'freq_E', 'freq_F', 'freq_G', 'freq_H', 'freq_I',
                       'freq_K', 'freq_L', 'freq_M', 'freq_N', 'freq_P', 'freq_Q',
                       'freq_R', 'freq_S', 'freq_T', 'freq_V', 'freq_W', 'freq_Y',
                       'ss_alpha', 'ss_beta', 'ss_coil',
                        'ss_turn',"GC_genomic",'kingdom',"phylum","name","count_protein"]

                gdf = gdf[cols].reset_index().dropna(subset=["kingdom"])
                gdf.to_csv(result_dir+"./df_" + ty + "_annotation_by_species.csv", index=False)


            gdf = pd.read_csv(result_dir+"./df_" + ty + "_annotation_by_species.csv")
            gdf.to_feather(result_dir+"./df_" + ty + "_annotation_by_species.feather")



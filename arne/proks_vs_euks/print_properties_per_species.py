

import os
import pandas as pd
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#props = ["iupred_long","iupred_short", "seg","top-idp","hessa","freq_S"]

props = ["iupred_long", "iupred_short", "top-idp","seg","hessa","length_translation"]
props += ["freq_A","freq_C","freq_D","freq_E","freq_F","freq_G","freq_H","freq_I","freq_K",
              "freq_L","freq_M","freq_N", "freq_P","freq_Q","freq_R","freq_S", "freq_T",
              "freq_V","freq_W","freq_Y"]



figures_dir = "./figures/"

dic_boxplot_limits = {"iupred_long":[0,50], "iupred_short":[0,50], 
                  "top-idp":[0,0.15],'length_translation':[0,800],
                  u'hessa':[0.85,1.2], 'seg':[0,25],"scampi":[0.02,0.1],
                "freq_S":[0.03,0.1],"freq_P":[0.01,0.08], 
                  "freq_I":[0.02,0.11],"freq_T":[0.03,0.08],
                          "freq_A":[0.00,0.18],
                          "freq_C":[0.0,0.06],
                          "freq_D":[0.02,0.1],
                          "freq_E":[0.01,0.12],
                          "freq_F":[0.00,0.1],
                          "freq_G":[0.03,0.15],
                          "freq_H":[0.00,0.06],
                          "freq_K":[0.00,0.12],
                          "freq_L":[0.05,0.15],
                          "freq_M":[0.01,0.05],
                          "freq_N":[0.00,0.1],
                          "freq_Q":[0.01,0.1],
                          "freq_R":[0.01,0.1],
                          "freq_V":[0.03,0.12],
                          "freq_W":[0.0,0.04],
                          "freq_Y":[0.00,0.07],
                          }



kingdoms = ["Eukaryota","Bacteria"]
datasets = ["full", "pfam", "domains", "linkers"]



dic_labels = {"GC":"GC%","GC_genomic":"GC% (genomic)",                'other':'Intermediate',                'ancient':'Ancient','orphan':'Orphans','group_orphan':"Genus\norphans",
    'hessa':"<Hydrophobicity> (Hessa)",                "top-idp":"TOP-IDP", "ss_alpha":"<Alpha propensity>","ss_beta":"<Beta propensity>","ss_coil":"<Coil propensity>",
    'ss_turn':"<Turn propensity>","length_translation":"Length (AA)",                'gc_third':'GC% (third nucleotide)','taxon_gc':'GC%        (taxonomic avg)','taxon_gc_r':'GC% (taxonomic avg)',
                'random':'Expected','iupred_long':'IUpred long(%AA)','length_aa':'<Length> (AA)','gc':'GC%',
                'iupred_short':'IUpred short(%AA)','iupred_cutoff':'IUpred cutoff (%AA)',
                'idp':'<TOP-IDP>','tm':'TM (%AA)','transmembrane':'TM(%AA)',"volatility":"<Volatility>",
                'iupred_short':'IUpred short(%AA)','iupred_cutoff':'IUpred cutoff(%AA)','seg':'<SEG>', 'hydropathy':'<Hydrophobicity>(Kyte-Doolittle)','gc':'GC%','aa_corr':'AA freq\ncorrelation with expected',
                 'aa_corr_dot':'AA freq\ncorrelation with expected','aa_corr_r':'AA freq\ncorrelation with expected', 'A': 'ALA', 'C':
                'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE',
                'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L':
                'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER',
                'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y':
                'TYR', 'freq_A': 'ALA', 'freq_C': 'CYS', 'freq_E': 'GLU',
                      'freq_D': 'ASP', 'freq_G': 'GLY', 'freq_F': 'PHE',
                      'freq_I': 'ILE', 'freq_H': 'HIS', 'freq_K': 'LYS',
                      'freq_M': 'MET', 'freq_L': 'LEU', 'freq_N': 'ASN',
                      'freq_Q': 'GLN', 'freq_P': 'PRO', 'freq_S':'SER', 'freq_R': 'ARG',
                      'freq_T': 'THR', 'freq_W':                'TRP', 'freq_V': 'VAL',
                      'freq_Y': 'TYR', 'scampi':'TM  (%AA)', 'slider':'Slider',
                      'disembl_coils':'DisEMBL                Coils (%AA)',
                      'disembl_rem465':"DisEMBL REM465 (%AA)",
                      'disembl_hotloops':"DisEMBL Hotloops (%AA)",
                      "globpipe":"GlobPlot (%AA)", "euk":"Eukaryota",
                "pro":"Bacteria",
                "arc":"Archaea",'Eukaryota':'Eukaryota','Bacteria':'Bacteria','Archaea':'Archaea'}




dfs = {}
for ty in datasets:
    dfs[ty] = pd.read_feather("./results/df_" + ty + "_annotation_by_species.feather")

phyla_euk = dfs["full"].loc[dfs["full"].kingdom == "Eukaryota"]["phylum"].dropna().unique()
phyla_prok = dfs["full"].loc[dfs["full"].kingdom == "Bacteria"]["phylum"].dropna().unique()
phyla_all = list(phyla_prok) + list(phyla_euk)


for prop in props:

    fig, ax = plt.subplots(figsize=(4,3))

    all_data = []
    labels = []
    for d in datasets:
        
        for k in kingdoms:

            kdf = dfs[d].loc[dfs[d].kingdom == k].dropna(subset=[prop])

            if prop != "length_translation":
                data = kdf[prop]/kdf["length_translation"]
            else:
                data = kdf[prop]

            all_data += [data]
            
            label = k + "_" + d
            labels += [label]


    plt.xticks(rotation=90)
    bp = ax.boxplot(all_data,labels=labels,patch_artist=True)
    
    
    for i,box in enumerate(bp['boxes']):
        if i%2 == 0:
            box.set(color='#93002F', linewidth=2)
            box.set(facecolor = '#E20048' )
        else:
            box.set(color='#00733E', linewidth=2)
            box.set(facecolor = '#00B060' )

    #~ pylab.setp(bp['boxes'],color='black',lw=2)
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['medians'], color='black')
    plt.setp(bp['fliers'], color='black', marker='+')
    
    ax.set_title(dic_labels[prop])
    #ax.set_ylim(dic_boxplot_limits[prop])
    

    figname = figures_dir + "boxplot_" + prop + ".png"
    print(figname)
    fig.savefig(figname,bbox_inches = "tight")
    plt.close(fig)



# phylum plots
for d in datasets:

    for prop in props:

        fig, ax = plt.subplots(figsize=(4,3))

        all_data = []
        labels = []
        
            
        for k in phyla_all:
          
            kdf = dfs[d].loc[dfs[d].phylum == k].dropna(subset=[prop])
 
            if prop != "length_translation":
                data = kdf[prop]/kdf["length_translation"]
            else:
                data = kdf[prop]

           
            if len(data) > 30:
                all_data += [data]
                
                label = k# + "_" + d
                labels += [label]


        plt.xticks(rotation=90)
        bp = ax.boxplot(all_data,labels=labels,patch_artist=True)
        
        for l,box in zip(labels,bp['boxes']):
            if l in set(phyla_euk):
                box.set(color='#93002F', linewidth=2)
                box.set(facecolor = '#E20048' )
            else:
                box.set(color='#00733E', linewidth=2)
                box.set(facecolor = '#00B060' )

        #~ pylab.setp(bp['boxes'],color='black',lw=2)
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['medians'], color='black')
        plt.setp(bp['fliers'], color='black', marker='+')


        ax.set_title(dic_labels[prop] + " - dataset: " + d)
        #ax.set_ylim(dic_boxplot_limits[prop])
        

        figname = figures_dir + "boxplot_phyla_" + prop + "_" + d + ".png"
        print(figname)
        fig.savefig(figname,bbox_inches = "tight")
        plt.close(fig)



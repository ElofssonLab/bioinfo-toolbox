#!/usr/bin/env python3
from __future__ import print_function

def get_aa_freq(aa,gc):
    codons = aa2codons[aa]

    prob_tot = []
    for codon in codons:
        prob_codon = []
        for N in codon:
            prob = gc/2
            if N in {"A","T","U"}:
                prob = 0.5-(gc/2)
            prob_codon += [prob]
            
        #print(codon, prob_codon,functools.reduce(mul,prob_codon))
        prob_tot += [functools.reduce(mul,prob_codon)] # Use functools.reduce() if you really need it; however, 99 percent of the time an explicit for loop is more readable.
    prob_tot = sum(prob_tot)
    
    return prob_tot



##########################################################################################################################

# for each amino acid, and separately for eukarya and bacteria, plot its frequency versus GC% for 4 curves:
# 1- average freq
# 2- a linear fit of the data
# 3- the frequency as it should be if only dependent on GC
# 4- a model function (w * prob_tot + k) with estimated parameters W and K
# note: for our big dataset (from NCBI) it is also possible to plot against the GC% of the CDS

# w and k are the parameters to optimize
def my_func(x,y,w,k):
    
    global aminoacid

    codons = aa2codons[aminoacid]
    
    res = []
    
    for gc in x:
        prob_tot = []
        for codon in codons:
            prob_codon = []
            for N in codon:
                prob = gc/2
                if N in {"A","T","U"}:
                    prob = 0.5-(gc/2)
                prob_codon += [prob]

            #print(codon, prob_codon,functools.reduce(mul,prob_codon))
            prob_tot += [functools.reduce(mul,prob_codon)]
        prob_tot = sum(prob_tot)

        res += [w * prob_tot + k]

    return res


def get_error(obs,exp):
    err = []
    for o,e in zip(obs,exp):
        err += [abs(o - e)]
    return np.mean(err)


def clamp(value, min_v, max_v):
    return max(min(value, max_v), min_v)

##########################################################################################################################

# define a function to plot a bar plot of the difference between AA frequencies in two datasets
# d1, d2 = the two dataset (could be the same dataset); datasets are "dataset","dataset-eggnog","dataset-pfam","dataset-pfam-full"
# k1, k2 = the kingdoms: can be either "Eukaryota" or "Bacteria"


def bar_plot_s1_vs_s2(d1,d2,k1,k2,uselog = True):
    fig, ax = plt.subplots(figsize=(6.,3))

    data_aa = []
    for aa in aas_sorted_topidp:
        prop = "freq_" + aa

        prop_values = {}
        prop_values_1 = dfs[d1].loc[dfs[d1].kingdom == k1][prop].mean()
        prop_values_2 = dfs[d2].loc[dfs[d2].kingdom == k2][prop].mean()

        if uselog == True:
            value = np.log(prop_values_1/prop_values_2)
        else:
            value = prop_values_1-prop_values_2
        data_aa += [value]

        
        
    ax.bar(range(len(data_aa)),data_aa, color="red")

    ylim = max([abs(v) for v in ax.get_ylim()])
    ylim = 0.03
    ax.set_ylim([-ylim, ylim])
    ax.set_xticks([0.0 + x for x in range(len(aas_sorted_topidp))])
    ax.set_xticklabels(aas_sorted_topidp)
    
    if uselog == True:
        ax.set_ylabel("log( freq(" + k1+ ") / freq(" + k2+ ") )")
    else:
        ax.set_ylabel("freq(" + k1+ ") - freq(" + k2+ ")")
    
    title = dic_labels[k1] + " (" + d1 + ")" + " vs "  + dic_labels[k2] + " (" + d2 + ")"
    ax.set_title(title)
    
    ax.yaxis.grid()
    ax.xaxis.grid()
    
    
    if uselog == True:
        figname = figures_dir + "aa_freq_change_log_odd_" + k1 + "_" + d1 + "_vs_" + k2 + "_" + d2 +  ".png"
    else:
        figname = figures_dir + "aa_freq_change_" + k1 + "_" + d1 + "_vs_" + k2 + "_" + d2 + ".png"
    
    print(figname)
    fig.savefig(figname,dpi=300)
#    fig.show()


def fmt_scientific(x):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def running_mean(x, N):
    cumsum = pd.np.cumsum(pd.np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N 



if __name__ == '__main__':


# coding: utf-8

# In[2]:

    import pandas as pd
    import os
    import functools
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats.stats import rankdata, square_of_sums, tiecorrect
    from scipy.stats import stats
    #from statsmodels.formula.api import ols
    from scipy.stats.stats import linregress
    #matplotlib.mathtext.SHRINK_FACTOR = 1
    
    #%matplotlib inline
    
    data_dir = "./data/"
    figures_dir = "./figures_new/"
    
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

# open the three reference files
    genomes_overview_file = data_dir + "overview.txt"
    genomes_euks_file = data_dir + "eukaryotes.txt"
    genomes_proks_file = data_dir + "prokaryotes.txt"
    
    df_overview = pd.read_csv(genomes_overview_file, sep="\t")

    df_proks = pd.read_csv(genomes_proks_file, sep = "\t")
    df_euks = pd.read_csv(genomes_euks_file, sep = "\t")
    df_reference = pd.concat([df_euks, df_proks])

    taxid2gc = df_reference.set_index("TaxID").to_dict()["GC%"]


    taxid2group = df_reference[["TaxID","Group"]].set_index("TaxID").to_dict()["Group"]
    group2kingdom = df_overview[["Group", "Kingdom"]].set_index("Group").to_dict()["Kingdom"]
    taxid2kingdom = {}

    for t in taxid2group:
        gr = taxid2group[t]
        taxid2kingdom[t] = group2kingdom[gr]
    
    taxid2name = df_reference[["#Organism/Name","TaxID"]].set_index("TaxID").to_dict()['#Organism/Name']
    name2taxid = df_reference[["#Organism/Name","TaxID"]].set_index("#Organism/Name").to_dict()['TaxID']


    # load random genes
    df_genes_random = pd.read_feather(data_dir + "df_genes_random.feather")


    # load datasets (per species)
    datasets = ["full", "pfam", "domains", "linkers"]
    #datasets += ["domains_exclusive_euk","domains_exclusive_bac","domains_orthologs_5"]

    # load datasets (per species)
    datasets_dom = ["full", "pfam", "domains", "linkers"]
    #datasets_prot = ["full","domains_exclusive_euk","domains_exclusive_bac","domains_orthologs_5"]
    
    dfs = {}
    for ty in datasets:   
        df = pd.read_feather("./results/df_" + ty+ "_annotation_by_species.feather")
        df.iupred_long = df.iupred_long*100
        df.iupred_short = df.iupred_short*100
        df.seg = df.seg*100
        dfs[ty] = df
        
    #~ dfs = {}
    #~ dfs["dataset"] = pd.read_feather(data_dir + "df_annotations_by_species.feather")
    #~ dfs["dataset-eggnog"] = pd.read_feather(data_dir + "df_annotations_eggnog_by_species.feather")
    #~ dfs["dataset-pfam"] = pd.read_feather(data_dir + "df_annotations_pfam_by_species.feather")
    #~ dfs["dataset-pfam-full"] = pd.read_feather(data_dir + "df_annotations_pfam_full_sequences_by_species.feather")


    # these taxa should be excluded because their genetic code is different, so translation is wrong
    mycotaxa = []
    for t in name2taxid:
        tax = t.split()[0]
        if tax in ["Mycoplasma","Spiroplasma","Ureaplasma","Mesoplasma"]:
            mycotaxa += [name2taxid[t]]
        

    # print some statistics
    for ds in datasets:
        print(ds)
        print(dfs[ds].kingdom.dropna().value_counts())
        print()
    


    # In[6]:

    # function to get a given AA frequency given the GC content
    aas_sorted_topidp = ['W','F','Y','I','M','L','V','N','C','T','A','G','R','D','H','Q','S','K','E','P']

    codons2aa = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        
    aa2codons = {}
    aa2codons["ALL"] = []
    for c in codons2aa:
        aa2codons[codons2aa[c]] = aa2codons.get(codons2aa[c],[]) + [c]
        if c != "_":
            aa2codons["ALL"] += [c]
            
    from operator import mul
    

    #get_aa_freq("ALL", 0.1)


    # In[7]:

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



    
    


    from scipy.optimize import curve_fit



    


    ##########################################################################################################################
    # for each amino acid, and separately for eukarya and bacteria, plot its frequency versus the taxonomic GC% for 4 curves:
    # 1- average freq
    # 2- a linear fit of the data
    # 3- the frequency as it should be if only dependent on GC
    # 4- a model function (w * prob_tot + k) with estimated parameters W and K
    # note: for our big dataset (from NCBI) it is also possible to plot against the GC% of the CDS

    # chose which dataset to use  for the following plots
    df = dfs["full"]
    
    #good_taxa = set(df.taxon_id.unique())
    good_taxa = set(df.taxon_id.unique()) - set(mycotaxa)

    res = {}

    for kingdom in ["Bacteria", "Eukaryota"]:
    
        res[kingdom] = {}
    
        for aa in [u'L', u'D', u'V', u'C', u'E', u'T', u'H', u'S',
            u'M', u'W', u'Y', u'F', u'Q', u'P', u'G', u'N', u'R', u'I',
            u'K', u'A']:
            
            global aminoacid
            aminoacid = aa
            
            # get the distribution of the frequency in bacteria
            d = df.loc[(df.kingdom == kingdom) & (df.taxon_id.isin(good_taxa))][["freq_"+aa,"GC_genomic"]].sort_values("GC_genomic").dropna(subset=["GC_genomic"])

            xdata = d["GC_genomic"]/100.0
            ydata = d["freq_"+aa]


            fig, ax = plt.subplots()
            
            # plot the data
            ax.scatter(xdata,ydata,color="lightgrey",edgecolor="None",alpha=0.3,label = "Data")
            
            
            # get an estimate of the error given only the average AA frequence
            ydata_mean = ydata.mean()
            ydata_std = ydata.std()
            avg_ydata = [ydata_mean for _ in ydata]
            ax.plot(xdata,avg_ydata,color="blue",linestyle="--",label = "Average")
            
            
            # create a first estimate of the fitting function and plot it
            initial_w = 1
            initial_k = 0
            estimated_ydata = my_func(xdata,ydata,initial_w,initial_k)
            ax.plot(xdata,estimated_ydata,color="darkred",linestyle="--",label = "GC Simple")
            
            
            # do a curve fit to estimate the parameters
            popt, pcov = curve_fit(my_func, xdata, ydata )
            _, est_w, est_k = popt
            # limit the parameters
            est_w = clamp(est_w, 0, 1)
            est_k = clamp(est_k, 0, ydata.mean())
            optimized_ydata = my_func(xdata, ydata, est_w, est_k)
            ax.plot(xdata, optimized_ydata, color="orange",lw=2, label='Function')
            
            
            # calculate and plot the fitting lines
            fit = np.polyfit(xdata,ydata,2)
            fit_fn = np.poly1d(fit) 
            fitted_ydata = fit_fn(xdata)
            plt.plot(xdata, fitted_ydata, lw=2,color="green",label="Linear Fit")
            
            
            res[kingdom][aa] = est_w
            
            
            # calculate and displays errors
            error_avg = get_error(ydata, avg_ydata) / ydata_std
            error_gc = get_error(ydata, estimated_ydata) / ydata_std
            error_func = get_error(ydata, optimized_ydata) / ydata_std
            error_fit = get_error(ydata, fitted_ydata) / ydata_std
            
            #text_str = "Avg Error = " + str(round(error_avg,4)) + "\n"
            #text_str += "GC Error = " + str(round(error_gc,4)) + "\n"
            #text_str += "Function Error = " + str(round(error_func/ydata_std,4)) + "\n"
            #text_str += "Fit Error = " + str(round(error_fit/ydata_std,4)) + "\n"
            text_str = "Error = " + str(round(get_error(fitted_ydata, optimized_ydata) / ydata_std,4)) + "\n"
            
            text_str += "W = " + str(round(est_w,5)) + "\n"
            text_str += "K = " + str(round(est_k,5)) 

            bb = ax.text(ax.get_xlim()[0], ax.get_ylim()[1], text_str, fontsize=12,
                    ha="left", va="top",
                    bbox=dict(boxstyle="square,pad=0",
                            ec=(0.0, 0.0, 0.0,1),
                            fc=(1., 1, 1,1),),    
                        )
            
            
            ax.set_xlabel("GC%")
            ax.set_ylabel("Frequency")
            ax.set_title(dic_labels[aa] +  " - " + dic_labels[kingdom])
            ax.legend(loc="lower left")
            
            
            fig_title =  figures_dir + "aafreq_" + aa + "_" + kingdom + ".eps"
            
            print(fig_title)
            fig.savefig(fig_title)
            #fig.show()
            
        #############################
        
    # examples

    
    for k in ["Bacteria", "Eukaryota"]:
            for j in ["Bacteria", "Eukaryota"]:
                for d in datasets:
                        for e in datasets:
                            bar_plot_s1_vs_s2(e,d, j, k,False)
                
#    bar_plot_s1_vs_s2("full","Pfam", "Eukaryota", "Eukaryota",False)
#    bar_plot_s1_vs_s2("full","Pfam", "Bacteria", "Bacteria",False)
#    bar_plot_s1_vs_s2("linkers","Pfam", "Eukaryota", "Eukaryota",False)
#    bar_plot_s1_vs_s2("linkers","Pfam", "Bacteria", "Bacteria",False)
#    bar_plot_s1_vs_s2("domains","Pfam", "Eukaryota", "Eukaryota",False)
#    bar_plot_s1_vs_s2("domains","Pfam", "Bacteria", "Bacteria",False)
#
#    bar_plot_s1_vs_s2("full", "full", "Eukaryota", "Bacteria",False)
#    bar_plot_s1_vs_s2("Pfam", "Pfam", "Eukaryota", "Bacteria",False)
#    bar_plot_s1_vs_s2("domains", "domains", "Eukaryota", "Bacteria",False)
#    bar_plot_s1_vs_s2("linkers", "linkers", "Eukaryota", "Bacteria",False)


    
    #################################################

    
    ######################################################

    # for a given dataset (ex. "dataset") plot two properties against each other

    df = dfs["full"]

    # the properties to plot
    properties_x = ["GC_genomic"]
    #properties_x = ["GC", "GC3", "GC_genomic"]

    properties_y = [u'top-idp', u'hessa', u'iupred_long', u'iupred_short', u'seg',
                    'ss_alpha', u'ss_beta', u'ss_coil', u'ss_turn',
        u'freq_A', u'freq_C', u'freq_D', u'freq_E', u'freq_F', u'freq_G',
        u'freq_H', u'freq_I', u'freq_K', u'freq_L', u'freq_M', u'freq_N',
        u'freq_P', u'freq_Q', u'freq_R', u'freq_S', u'freq_T', u'freq_V',
        u'freq_W', u'freq_Y']
        
        
        
        
    dic_limits={"hessa":[0.85,1.25],"GC":[20,80],"GC_genomic":[20,80],"taxon_gc":[30,65],"taxon_gc_r":[30,65],"gc_third":[30,65],
                'aa_corr':[0.8,1],'aa_corr_dot':[0.8,1.0],'aa_corr_r':[0.4,1],
                "top-idp":[0,0.15], 'W':[0.006,0.02], 'C':[0.01,0.04],
                'H':[0.01,0.04], 'M':[0.01,0.04], 'Y':[0.01,0.06],
                'Q':[0.02,0.05], 'P':[0.02,0.11], 'F':[0.01,0.08],
                'R':[0.02,0.1], 'G':[0,0.12], 'D':[0.02,0.06],
                'A':[0,0.15], 'V':[0.04,0.08], 'N':[0.01,0.08],
                'T':[0.04,0.08], 'E':[0.02,0.08], 'I':[0,0.12],
                'K':[0.01,0.1], 'S':[0.05,0.11], 'L':[0.06,0.12],
                'freq_W':[0.0,0.025], 'freq_C':[0.00,0.04],
                'freq_H':[0.01,0.035], 'freq_M':[0.01,0.065],
                'freq_Y':[0.01,0.06], 'freq_Q':[0.01,0.07],
                'freq_P':[0.0,0.08], 'freq_F':[0.018,0.08],
                'freq_R':[0.02,0.11], 'freq_G':[0.02,0.12],
                'freq_D':[0.02,0.1], 'freq_A':[0,0.18],
                'freq_V':[0.04,0.11], 'freq_N':[0.01,0.095],
                'freq_T':[0.03,0.08], 'freq_E':[0.02,0.1],
                'freq_I':[0.02,0.14], 'freq_K':[0.01,0.13],
                'freq_S':[0.04,0.11], 'freq_L':[0.075,0.15],
                "iupred_long":[0,40], "tm":[0,12],'idp':[0,0.20],
                "seg":[0,25],"length_aa":[0,720],'iupred_short':[0,40],
                'iupred_cutoff':[0,50],'volatility':[0.7,0.85],'hydropathy':[-1,1],
                "ss_alpha":[-0.05,0.03],"ss_beta":[-0.1,0.02],"ss_coil":[-0.03,0.02],"ss_turn":[-0.11,0.0],
                'slider':[10,80],'disembl_coils':[40,100],'disembl_rem465':[0,80],'disembl_hotloops':[10,90],
                "globpipe":[0,100] } # hardcoded axis limits


    dic_colors = {"Eukaryota":"#9daf6a", "Bacteria":"#b45a5c", "Archaea":"#377eb8"}
    dic_colors_fit = {"Eukaryota":"#bdcf8a", "Bacteria":"#f43a3c", "Archaea":"#579ed8"}
    dic_markers = {"Eukaryota":"o", "Bacteria":"x", "Archaea":"o"}
        

    for prop_1 in properties_x:
        for prop_2 in properties_y:

            fig, ax = plt.subplots(figsize=(6,6))
        
            # plot the values as a scatter plot
            prop_1_values = {}
            prop_2_values = {}
            for kingdom in ['Eukaryota','Bacteria']:
                

                good_taxa = set(df.taxon_id.unique()) - set(mycotaxa)
                
                temp_df = df.loc[(df.kingdom == kingdom) & (df.taxon_id.isin(good_taxa))][list(set([prop_1,prop_2,"iupred_short"]))].dropna()
                
                prop_1_values[kingdom] = temp_df[prop_1].values
                prop_2_values[kingdom] = temp_df[prop_2].values
                
                ax.scatter(prop_1_values[kingdom],prop_2_values[kingdom], edgecolor="none", facecolor=dic_colors[kingdom],alpha = 1,s=50,label=dic_labels[kingdom])
                #sct = ax.scatter(prop_1_values[kingdom],prop_2_values[kingdom],marker=dic_markers[kingdom], c=colors_list, edgecolor="none",alpha = 1,s=70,label=dic_labels[kingdom])
                
                # calculate and plot the fitting lines
                fit = np.polyfit(list(prop_1_values[kingdom]),list(prop_2_values[kingdom]),1)
                fit_fn = np.poly1d(fit) 
                plt.plot(prop_1_values[kingdom], fit_fn(prop_1_values[kingdom]), '--k',lw=2,color=dic_colors_fit[kingdom])

                #slope, intercept, r_value, p_value, std_err = linregress(prop_1_values, prop_2_values)
                #print(slope, intercept, r_value, p_value, std_err)
                

            # plot the random running average
            if prop_1 in df_genes_random.columns and prop_2 in df_genes_random.columns:
                prop_1_random = running_mean(pd.np.array(df_genes_random.sort_values(prop_1)[prop_1].astype(float)),500)
                prop_2_random = running_mean(pd.np.array(df_genes_random.sort_values(prop_1)[prop_2].astype(float)),500)
                ax.plot(prop_1_random,prop_2_random,color="black",label="Random")
                
            # plot the theoretical line
            if prop_2.split("_")[0] == "freq" and prop_1 in ["GC","GC3","GC_genomic"]:
                aa = prop_2.split("_")[1]
                p = []
                for gc in np.arange(0.1,0.9,0.001):
                    p += [get_aa_freq(aa, gc)]
                
                ax.plot(np.arange(10,90,0.1),p,color="blue",lw=2,label="Theoretical")


            # set the limits
            if prop_1 in dic_limits:
                ax.set_xlim(dic_limits[prop_1])
                
            if prop_2 in dic_limits:
                ax.set_ylim(dic_limits[prop_2])
                
                
            # set the labels
            ax.set_xlabel(dic_labels.get(prop_1,prop_1))
            ax.set_ylabel(dic_labels.get(prop_2,prop_2))
            
            # plot legend
            ax.legend(loc="upper right")
            
            # ranksums tests
            ranksum_string = ""
            ranksum_string += "Eukarya vs Bacteria = " + fmt_scientific(stats.ranksums(prop_2_values["Eukaryota"], prop_2_values["Bacteria"])[1]) + "\n"
            #ranksum_string += "Eukarya vs Archaea = " + fmt_scientific(stats.ranksums(prop_2_values["Eukarya"], prop_2_values["Archaea"])[1]) + "\n"
            #ranksum_string += "Bacteria vs Archaea = " + fmt_scientific(stats.ranksums(prop_2_values["Bacteria"], prop_2_values["Archaea"])[1]) + "\n" + "\n"
            
            # Pearson correlation coefficient
            ranksum_string += "Eukarya R = " + str(round(stats.pearsonr(prop_1_values["Eukaryota"], prop_2_values["Eukaryota"])[0],2)) + "\n"
            ranksum_string += "Bacteria R = " + str(round(stats.pearsonr(prop_1_values["Bacteria"], prop_2_values["Bacteria"])[0],2)) + "\n"
            #ranksum_string += "Archaea R = " + str(round(stats.pearsonr(prop_1_values["Archaea"], prop_2_values["Archaea"])[0],2))
            
            
            # calculate ABS difference between Euk-theoretical and Prok-theoretical
            if prop_2.split("_")[0] == "freq" and prop_1 in ["GC","GC3","GC_genomic"]:
                
                diff_from_theo = {}
                
                aa = prop_2.split("_")[1]
                for kingdom in ['Eukaryota','Bacteria']:
                    p = []
                    for vx,vy in zip (prop_1_values[kingdom], prop_2_values[kingdom]):
                        p += [abs(vy-get_aa_freq(aa, vx/100.0))]
                        
                    ranksum_string += "abs(" + kingdom + " - theoretical) = " + str(round(np.mean(p),4))+ "\n"
                    
                    
            bb = ax.text(ax.get_xlim()[0], ax.get_ylim()[1], ranksum_string, fontsize=12,
                ha="left", va="top",
                bbox=dict(boxstyle="square,pad=0",
                        ec=(0.0, 0.0, 0.0,1),
                        fc=(1., 1, 1,1),),    
                    )
            
        
            
            figname = figures_dir + prop_2 + "_vs_" + prop_1 + "_pfam.eps"
            print(figname)
            fig.savefig(figname)
            #fig.show()
            
            ##########################################################################################################################
            
            
            
            
        ##########################################################################################################################
        
    # plot a series of properties as boxplots, with each box corresplonding to a dataset-kingdom combination
    
    # the properties to plot
    props = ["iupred_long", "iupred_short", "top-idp","seg","length_translation"]
    props += ["freq_A","freq_C","freq_D","freq_E","freq_F","freq_G","freq_H","freq_I","freq_K",
                  "freq_L","freq_M","freq_N", "freq_P","freq_Q","freq_R","freq_S", "freq_T",
                  "freq_V","freq_W","freq_Y"]



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


    for prop in props:
    
        all_data = []
        labels = []
        for d in datasets:
            for k in kingdoms:

                data = dfs[d].loc[dfs[d].kingdom == k][prop].dropna()

                all_data += [data]
                
                label = k + "_" + d
                labels += [label]
                
        fig, ax = plt.subplots(figsize=(4,3))
        
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
        ax.set_ylim(dic_boxplot_limits[prop])
        
        
        figname = figures_dir + "boxplot_" + prop + ".png"
        print(figname)
        fig.savefig(figname,bbox_inches = "tight")
        #fig.show()

         
        #************************************************************************
        # Not we divide this into subfigs
        all_data = []
        labels = []
        for d in datasets_dom:
            for k in kingdoms:

                data = dfs[d].loc[dfs[d].kingdom == k][prop].dropna()
                if len(data) > 1:
                    all_data += [data]
					
                    label = k + "_" + d
                    labels += [label]
					
        fig, ax = plt.subplots(figsize=(4,3))
        
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
        ax.set_ylim(dic_boxplot_limits[prop])
        
        
        figname = figures_dir + "boxplot2_" + prop + ".png"
        print(figname)
        fig.savefig(figname,bbox_inches = "tight")        
        ############################################



        '''
        #************************************************************************
        # Not we divide this into subfigs
        all_data = []
        labels = []
        for d in datasets_prot:
            for k in kingdoms:

                data = dfs[d].loc[dfs[d].kingdom == k][prop].dropna()
                if len(data) > 1:
                    all_data += [data]
					
                    label = k + "_" + d
                    labels += [label]
                
        fig, ax = plt.subplots(figsize=(4,3))
        
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
        ax.set_ylim(dic_boxplot_limits[prop])
        
        
        figname = figures_dir + "boxplot3_" + prop + ".png"
        print(figname)
        fig.savefig(figname,bbox_inches = "tight")        
        ############################################
        '''



    '''
    # plot a series of properties as boxplots, with each box corresplonding to phylum
    
    # the properties to plot
    props = ["iupred_long", "iupred_short", "top-idp",'length_translation', u'hessa', 'seg','GC_genomic']
    props += ["freq_A","freq_C","freq_D","freq_E","freq_F","freq_G","freq_H","freq_I","freq_K",
                  "freq_L","freq_M","freq_N", "freq_P","freq_Q","freq_R","freq_S", "freq_T",
                  "freq_V","freq_W","freq_Y"]



    phyla_euk = dfs["full"].loc[dfs["full"].kingdom == "Eukaryota"]["phylum"].unique()
    phyla_proks = dfs["full"].loc[dfs["full"].kingdom == "Bacteria"]["phylum"].unique()
    all_phyla = list(phyla_proks) + list(phyla_euk)
    
    
    
    for prop in props:    
        data = []
        labels = []
        for phylum in all_phyla:
            dd = dfs["full"].loc[dfs["full"].phylum == phylum][prop].dropna()
            if len(dd) > 10:
                data += [dd]
                labels += [phylum]
        
        fig, ax = plt.subplots(figsize=(7,4))

        plt.xticks(rotation=90)
        bp = ax.boxplot(data,labels=labels,patch_artist=True)
        
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
        
    
        ax.set_title(dic_labels[prop])
        
        
        figname = figures_dir + "boxplots_by_phylum_" + prop + ".png"
        print(figname)
        fig.savefig(figname,bbox_inches = "tight")
        #fig.show()
        
        ##################################################
    '''





    #~ # plot contribution of
    #~ props = ["iupred_long", "iupred_short"]
    #~ l = "length_translation"
    #~ all_data = []
    #~ for k in kingdoms:
        #~ for d in datasets:
            #~ length = dfs[d].loc[dfs[d].kingdom == k][l].dropna()
            #~ print(k,d,"length",length.mean())

            #~ for p in props:
                #~ data = dfs[d].loc[dfs[d].kingdom == k][p].dropna()
                #~ length = dfs[d].loc[dfs[d].kingdom == k][l].dropna()
                #~ weight = data*length
                #~ print(k,d,p,weight.mean()/100,weight.mean()/(100*length.mean()))

    
    
    # plot the top-idp score of each AA weighted by the average frequency
    # of that AA in each kingdom, for each dataset
    
    # Top-IDP scale
    top_idp = {'A':0.06, 'C' :  0.02, 'D' : 0.192, 'E' : 0.736,
        'F' :  -0.697, 'G' : 0.166, 'H':0.303, 'I' :  -0.486,
        'K' : 0.586, 'L' :  -0.326, 'M': -0.397, 'N' : 0.007,
        'P' : 0.987, 'Q' : 0.318, 'R' : 0.180, 'S':  0.341,
        'T' : 0.059, 'V' :  -0.121, 'W':  -0.884, 'Y' : -0.510}

    for d in datasets:
        weight_euk = {}
        weight_bac = {}
        weight_diff = {}
        
        
        for aa in aas_sorted_topidp:
            if d not in (["domains_exclusive_euk","domains_exclusive_bac"]):
                weight_euk[aa] = dfs[d].loc[dfs[d].kingdom == "Eukaryota"]["freq_"+aa].mean() * top_idp[aa]
                weight_bac[aa] = dfs[d].loc[dfs[d].kingdom == "Bacteria"]["freq_"+aa].mean() * top_idp[aa]
            else:
                weight_euk[aa] = dfs["domains_exclusive_euk"]["freq_"+aa].mean() * top_idp[aa]
                weight_bac[aa] = dfs["domains_exclusive_bac"]["freq_"+aa].mean() * top_idp[aa]
                
            weight_diff[aa] = weight_euk[aa] - weight_bac[aa]

            
        fig, ax = plt.subplots()

        #ax.plot(range(len(aas_sorted_topidp)),[weight_euk[aa] for aa in aas_sorted_topidp],color = "red",marker="o",label="Euk") 
        #ax.plot(range(len(aas_sorted_topidp)),[weight_bac[aa] for aa in aas_sorted_topidp],color = "blue",marker="^",label="Bac") 
        ax.bar(range(len(aas_sorted_topidp)),[weight_diff[aa] for aa in aas_sorted_topidp],color = "green",label="Diff") 
        
        ax.set_xticklabels(aas_sorted_topidp)
        ax.set_xticks([p for p in range(len(aas_sorted_topidp))])

        ax.set_ylabel("TOP-IDP * frequency\n(Euk - Bac)")
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)

        #ax.legend(loc="upper left")
        ax.set_ylim([-0.006, 0.01])
        
        if d not in (["domains_exclusive_euk","domains_exclusive_bac"]):
            ax.set_title("Dataset: " + d)
            fname = figures_dir + "top-idp-weight_diff_" + d + ".eps"
        else:
            ax.set_title("Dataset: domains_exclusive")
            fname = figures_dir + "top-idp-weight_diff_domains_exclusive.eps"
        
        print(fname)
        fig.savefig(fname,bbox_inches = "tight")
        #plt.show()

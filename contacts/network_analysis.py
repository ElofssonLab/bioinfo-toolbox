import sys
import math
import os

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pylab import cm

import parse_contacts
import parse_fasta
import parse_mcl


def clustering(cfile_name, reflen):
    gfile_name = '%s.abc' % '.'.join(cfile_name.split('.')[:-1])
    mclfile_name = '%s.mcl' % '.'.join(cfile_name.split('.')[:-1])
    cfile = open(cfile_name, 'r')
    clist = parse_contacts.parse(cfile)
    G = nx.Graph()    

    """
    # residue based graph
    for i, c in enumerate(clist):
        aa1 = int(c[1])
        aa2 = int(c[2])
        dist = abs(aa1 - aa2)
        G.add_edge(aa1, aa2, weight=dist)
        if i > reflen:
            break
    """
    # get top L contacts
    clist_top = []
    for i, c in enumerate(clist):
        aa1 = int(c[1])
        aa2 = int(c[2])
        dist = abs(aa1 - aa2)
        if dist < 12:
            continue
        clist_top.append(c)
        if i >= reflen:
            break

    # create graph
    gfile = open(gfile_name, 'w')
    for c1 in clist_top:
        c1_name = '%d-%d' % (c1[1], c1[2])
        c1_pos = (c1[1], c1[2])
        for c2 in clist_top:
            if c1 == c2:
                continue
            c2_name = '%d-%d' % (c2[1], c2[2])
            c2_pos = (c2[1], c2[2])
            dist = np.linalg.norm(np.array(c2_pos) - np.array(c1_pos))

            # naive clustering (distance cutoff)
            # networkx object G = naive clustered graph
            if dist > 20:
                G.add_node(c1_pos)
                G.add_node(c2_pos)
            else:
                gfile.write('%s\t%s\t%s\n' % (c1_pos, c2_pos, dist))
                G.add_edge(c1_pos, c2_pos, weight=dist)
            #G.add_edge(c1_pos, c2_pos, weight=dist)
            
            # write file for mcl clustering
            #gfile.write('%s\t%s\t%s\n' % (c1_pos, c2_pos, dist))
            #gfile.write('%s\t%s\n' % (c1_pos, c2_pos))
    gfile.close()

    # run mcl clustering on graph file
    os.system('mcl %s --abc -I 4.5 --analyze=y -V all -o %s' % (gfile_name, mclfile_name))

    # read mcl output
    G_mcl = parse_mcl.get_clustering(open(mclfile_name, 'r'), G)
    mcl_stats = parse_mcl.get_stats(open(mclfile_name, 'r'))

    #print G_mcl.number_of_nodes()
    #print G_mcl.number_of_edges()
    #print G_mcl.nodes()
    #print G_mcl.edges()

    #print mcl_stats
    stat_names = ['clusters', 'max', 'ctr', 'avg', 'min', 'DGI', 'TWI', 'TWL', 'sgl', 'qrt', 'efficiency', 'massfrac', 'areafrac']
    stat_str = ''
    for stat in stat_names:
        stat_str += '%s\t' % mcl_stats[stat]
    print stat_str

    #print G.edges(data=True)
    clust_coeff = nx.average_clustering(G_mcl, weight='weight')
    #print clust_coeff
    #print len(clust_coeff)
    min_span_tree = nx.minimum_spanning_tree(G_mcl, weight='weight')
    #pos = dict(zip(G,G))
    #nx.draw(G, pos, node_size=10, with_labels=False)
    pos_mcl = dict(zip(G_mcl,G_mcl))
    #nx.draw(G_mcl, pos_mcl, node_size=10, with_labels=False)
    #plt.savefig('%s.pdf' % '_'.join(cfile_name.split('/')[-2:]))
    #plt.show()
    #return mcl_stats
    return G_mcl



if __name__ == '__main__':

    #cfile_name = sys.argv[1]
    #reflen = int(sys.argv[2])

    stat_names = ['clusters', 'max', 'ctr', 'avg', 'min', 'DGI', 'TWI', 'TWL', 'sgl', 'qrt', 'efficiency', 'massfrac', 'areafrac']
    stat_dict = {}

    idfile = open(sys.argv[1], 'r')
    outfile = open(sys.argv[2], 'w')
    outfile.write('ID\t%s\n' % '\t'.join(stat_names))
    for line in idfile:
        id = line.strip().split()
        if len(id) == 1:
            id = id[0]
            seqfile = open('data/%s:A/sequence.fa' % id, 'r')
            seq = parse_fasta.read_fasta(seqfile).values()[0][0]
            reflen = len(seq)
            #cfile_name = 'data/%s:A/pconse.results' % id
            cfile_name = 'data/%s:A/mar23-hh0hh4hh10hhb40jh0jh4jh10jhm40-psicovplmdca.results' % id
            #stat_lst.append(clustering(cfile_name, reflen))
            stat_dict[id] = clustering(cfile_name, reflen)
            stat_str = '%s\t' % id
            for stat in stat_names:
                stat_str += '%s\t' % stat_dict[id][stat]
            stat_str = stat_str.strip()
            stat_str += '\n'
            outfile.write(stat_str)
        else:
            seqfile = open('pfam/data/pconsc_predictions/%s/%s.fa' % (id[0], id[1]), 'r')
            seq = parse_fasta.read_fasta(seqfile).values()[0][0]
            reflen = len(seq)
            cfile_name = 'pfam/data/pconsc_predictions/%s/%s.out' % (id[0], id[1])
            stat_dict[id[0]] = clustering(cfile_name, reflen)
            stat_str = '%s\t' % id[0]
            for stat in stat_names:
                stat_str += '%s\t' % stat_dict[id[0]][stat]
            stat_str = stat_str.strip()
            stat_str += '\n'
            outfile.write(stat_str)

    idfile.close()
    outfile.close()
    """    
    outfile = open(sys.argv[2], 'w')
    outfile.write('\t'.join(stat_names) + '\n')
    for stat_dict in stat_dict:
        stat_str = ''
        for stat in stat_names:
            stat_str += '%s\t' % stat_dict[stat]
        stat_str = stat_str.strip()
        stat_str += '\n'
        outfile.write(stat_str)
    """
    outfile.close()



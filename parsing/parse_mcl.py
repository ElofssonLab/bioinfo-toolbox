import sys

import networkx as nx


def get_curr_clust(node_lst, G):

    curr_clust = nx.Graph()

    for node1 in node_lst:
        for node2 in node_lst:
            #if node1 == node2:
            #    continue
            if (node1, node2) in G.edges():
                curr_clust.add_edge(node1, node2, G[node1][node2])
                #print curr_clust.edges()
    
    return curr_clust


def get_clustering(infile, G):

    G_mcl = nx.Graph()

    for line in infile:
        if len(line.split('=')) == 1:
            line_arr = line.split('\t')
            if len(line_arr) <= 1:
                continue
            node_lst = map(eval, line_arr)
            curr_clust = get_curr_clust(node_lst, G)
            G_mcl.add_edges_from(curr_clust.edges(data=True))
    infile.close()

    return G_mcl


def get_stats(infile):

    mcl_stats = {}

    for line in infile:
        if len(line.split('=')) > 1:
            line_arr = line.split()
            for stat in line_arr:
                name = stat.split('=')[0]
                val = stat.split('=')[1]
                mcl_stats[name] = val

    infile.close()

    return mcl_stats

#!/usr/bin/env python

import sys
import argparse

try:
    import graphviz as gv
except:
    print "The graphviz library needs to be installed"
    sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="the hmg filename")
parser.add_argument("out", help="image output name")
parser.add_argument("-if", "--imageformat", help="output image format",
                    default="png")
parser.add_argument("-we", "--witoutclusteredges", help="draw without edges" +
                    " between clusters", action="store_true")
parser.add_argument("-debug", help="Display image automatically",
                    action="store_true")

args = parser.parse_args()
filename = args.filename
outname = args.out
imageFormat = args.imageformat
view = True if args.debug else False

try:
    hmgFile = open(filename)
    fullFile = hmgFile.read()
    hmgFile.close()
except:
    print "File " + filename + " was not found"
    sys.exit(1)
fileList = fullFile.split('\n\n\n')  # 0 = header, 1 = body (all modules)


headerList = fileList[0].split('\n')

nameRow = 0
while True:
        if headerList[nameRow].startswith('NAME'):
            break
        nameRow += 1


bodyList = fileList[1].split('---------------' +
                             '----------------------------------------')

colorDict = {"Singlenode": "blue",
             "Highway": "grey",
             "U_Turn": "tomato",
             "Singleloop": "darkorange"}
i = 0
graphs = {}
headGraph = gv.Digraph(engine="dot", name="headGrapher", format=imageFormat,
                       graph_attr={},
                       node_attr={},
                       edge_attr={})

edges = []
firstNodes = []
for module in bodyList[:-1]:

    cluster = gv.Digraph(name="cluster"+str(i),
                         node_attr={},
                         graph_attr={'fontsize': '24'})
    i += 1
    nodes = []
    targets = []
    first = True
    for line in module.split('\n'):
        if line.startswith('Module'):
            module = line.split(': ')[1]
        elif line.startswith('Type'):
            modType = line.split(': ')[1]
            cluster.attr('graph', _attributes={"label": modType+'\n'+module})
        elif line.startswith('Vertex type'):
            vertexType = line.split(': ')[1]
        elif line.startswith('Vertex label'):
            vertexLabel = line.split(': ')[1]
        elif line.startswith('Vertex'):
            source = line.split(' ')[1][:-1]
            nodes.append(source)
            if module == 's':
                color = "green"
            elif module == 'e':
                color = "red"
            else:
                color = colorDict[modType]
            cluster.node(source, label='Vertex ' + source, _attributes={
                'style': 'filled',
                'fillcolor': color,
                'labelfontsize': '16'})
        elif line.startswith('Transition prob'):
            transition = 'transition'
        elif line.startswith('End transition'):
            transition = 'end_transition'
        elif line.startswith('\tVertex'):
            target = line.strip().split(': ')[0].split(' ')[1]
            if module == 's':
                color = "green"
            elif module == 'e':
                color = "red"
            else:
                color = colorDict[modType]

            targets.append((source, target, color))

    for source, target, color in targets:
        if target in nodes:
            cluster.edge(source, target, _attributes={'color': color})
        else:
            edges.append((source, target, color))

    headGraph.subgraph(cluster)

if not args.witoutclusteredges:
    for edge in edges:
        source, target, color = edge
        headGraph.edge(source, target, _attributes={'constraint': 'false',
                                                    'color': color})
headGraph.render(filename=outname, view=view, cleanup=True)

import sys
from xml.etree.ElementTree import Element, SubElement
from xml.etree import ElementTree
from xml.dom import minidom


# Pretty xml printing
def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


# Filename from argument, read full file and split header from body
filename = sys.argv[1]

hmgFile = open(filename)
fullFile = hmgFile.read()
hmgFile.close()
fileList = fullFile.split('\n\n\n')  # 0 = header, 1 = body (all modules)

# Start xml document
graphml = Element('graphml',
                  attrib={"xmlns": "http://graphml.graphdrawing.org/xmlns",
                          "xmlns:xsi":
                          "http://www.w3.org/2001/XMLSchema-instance",
                          "xsi:schemaLocation":
                          "http://graphml.graphdrawing.org/xmlns " +
                          "http://graphmlgraphdrawing.org/xmlns/" +
                          "1.0/graphml.xsd"})

headerList = fileList[0].split('\n')

nameRow = 0
while True:
        if headerList[nameRow].startswith('NAME'):
            break
        nameRow += 1

graphName = headerList[nameRow].split(': ')[1]

graph = SubElement(graphml, 'graph', attrib={"id": graphName,
                                             "edgedefault": "directed"})
numberOfVerticesRow = nameRow
while True:
    if headerList[numberOfVerticesRow].startswith('NR OF VERTICES'):
        break
    numberOfVerticesRow += 1

numberOfVertices = headerList[numberOfVerticesRow].split(': ')[1]

for i in range(int(numberOfVertices)):
    node = SubElement(graph, 'node', attrib={"id": str(i)})

bodyList = fileList[1].split('\n')

transition = 'starting'
for line in bodyList:
    if line.startswith('Module'):
        module = line.split(': ')[1]
    if line.startswith('Type'):
        modType = line.split(': ')[1]
    if (line.startswith('Vertex') and
            not line.startswith('Vertex type') and
            not line.startswith('Vertex label')):
        source = line.split(' ')[1][:-1]
    if line.startswith('Vertex type'):
        vertexType = line.split(': ')[1]
    if line.startswith('Vertex label'):
        vertexLabel = line.split(': ')[1]

    if line.startswith('\tVertex'):
        target = line.strip().split(': ')[0].split(' ')[1]
        if transition == 'transition':
            edge = SubElement(graph, 'edge', attrib={"source": source,
                                                     "target": target})
    if line.startswith('Transition prob'):
        transition = 'transition'
    if line.startswith('End transition'):
        transition = 'end_transition'

# print prettify(graphml)

outFile = open(filename + '.xml', 'w')
outFile.write(prettify(graphml))
outFile.close()

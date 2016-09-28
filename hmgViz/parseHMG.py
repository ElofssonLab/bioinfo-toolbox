import sys
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.etree import ElementTree
from xml.dom import minidom

### Pretty xml printing
def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


## Filename from argument, read full file and split header from body
filename = sys.argv[1]

hmgFile = open(filename)
fullFile = hmgFile.read()
hmgFile.close()
fileList = fullFile.split('\n\n\n') # 0 = header, 1 = body (all modules)

## Start xml document
root = Element('root')
root.append(Comment('Generated from ' + filename))

## Generate up the first child, the header
head = SubElement(root, 'head')
header = fileList[0]

for row in header.split('\n')[1:]:
    keyvalues = row.strip().split(': ')
    if len(keyvalues)>1:
        key, value = keyvalues
        key = key.replace(' ','_')
        key = SubElement(head, key)
        key.text = value
    else:
        key = keyvalues[0].replace(' ','_')[:-1]
        key = SubElement(head, key)


## Prepare the body, take away the modules header line
i = fileList[1].index('\n')
modules = fileList[1][i+1:]

## Split individual modules into list
moduleList = modules.split('-------------------------------------------------------')[:-1]

body = SubElement(root, 'body')

## Run through each module block
for moduleBlock in moduleList:
    module = SubElement(body, 'module')
    subModules = moduleBlock.strip().split('\n\n')

    ## First we have the modules header
    for row in subModules[0].split('\n'):
        keyvalues = row.strip().split(': ')
        if len(keyvalues)>1:
            key, value = keyvalues
            key = key.replace(' ','_')
            key = SubElement(module, key)
            key.text = value
        else:
            key = keyvalues[0].replace(' ','_')[:-1]
            key = SubElement(module, key)
    vertex = SubElement(module, 'vertex')
    ## For each vertex...
    for vertexBlock in subModules[1:]:
        ## First the vertex name
        vertexLines = vertexBlock.split('\n')
        name = SubElement(vertex, 'name')
        name.text = vertexLines[0].split(' ')[1][:-1]

        ## Then the first part of the rest of the header
        for row in vertexLines[1:5]:
            keyvalues = row.strip().split(': ')
            if len(keyvalues)>1:
                key, value = keyvalues
                key = key.replace(' ','_')
                key = SubElement(vertex, key)
                key.text = value
            else:
                key = keyvalues[0].replace(' ','_')[:-1]
                key = SubElement(vertex, key)

        ## Second part, '=' instead of ':'
        for row in vertexLines[5:8]:
            keyvalues = row.strip().split(' = ')
            if len(keyvalues)>1:
                key, value = keyvalues
                key = key.replace(' ','_')
                key = SubElement(vertex, key)
                key.text = value
            else:
                key = keyvalues[0].replace(' ','_')[:-1]
                key = SubElement(vertex, key)

        ## Transitions
        transitions = SubElement(vertex, 'transitions')
        i = 9 ## Startline withing the vertexline
        while True:
            if vertexLines[i].startswith(' ') or vertexLines[i].startswith('\t'):
                #print 'right here'
                vertexNr, probs = vertexLines[i].split(': ')
                nr = SubElement(transitions, 'V' + str(vertexNr.strip().split(' ')[1]))
                nr.text = probs
            else:
                break
            i += 1

        ## End Transitions
        endtransitions = SubElement(vertex, 'endtransitions')
        i += 1 ## SOne step forward to start with the correct line
        while True:
            if vertexLines[i].startswith(' ') or vertexLines[i].startswith('\t'):
                #print 'right here'
                vertexNr, probs = vertexLines[i].split(': ')
                nr = SubElement(endtransitions, 'V' + str(vertexNr.strip().split(' ')[1]))
                nr.text = probs
            else:
                break
            i += 1

        ## Emission
        emissions = SubElement(vertex, 'emissions')
        i += 1 ## SOne step forward to start with the correct line
        while True:
            if vertexLines[i].startswith(' ') or vertexLines[i].startswith('\t'):
                #print 'right here'
                aa, probs = vertexLines[i].split(': ')
                nr = SubElement(emissions, aa.strip())
                nr.text = probs
            else:
                break
            i += 1
            if i>=len(vertexLines):
                break

    #print prettify(root)
    #break


print prettify(root)
outFile = open(filename + '.out', 'w')
outFile.write(prettify(root))
outFile.close()

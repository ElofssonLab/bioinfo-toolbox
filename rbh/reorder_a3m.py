import numpy as np
import sys

'''
To trim the two interacting proteins' output multiple sequence alignment of jackhmmer; as some sequences might lose under certain
Evalue and iteration settings. So this script trim both the sequences(ordered&same length fasta file) and MSA(a3m files[one liner]) 
to make sure that sequences and MSA are matching. 

input:
 pr1/pr2: fasta sequences of interacting protein 1/2.
 msa1/msa2: a3m file(one liner) of the protein 1/2.
 
output:
   
'''
pr1 = sys.argv[1]
pr2 = sys.argv[2]
msa1= sys.argv[3]
msa2 = sys.argv[4]

with open(pr1,'r') as f:
    fasta1=f.readlines()
    fasta1=[x.strip() for x in fasta1]
    order1 = []
    for line in fasta1:
        if line.startswith(">"):
            order1.append(line.split()[0][1:])

with open(pr2,'r') as f:
    fasta2=f.readlines()
    fasta2=[x.strip() for x in fasta2]
    order2 = []
    for line in fasta2:
        if line.startswith(">"):
            order2.append(line.split()[0][1:])

#warning: two fasta sequences must be ordered and in same length!

if len(order1)!=len(order2):
    print('fatal error: two fasta files do not match!')
    sys.exit("errors!")

with open(msa1, 'r') as a3m:
    a3m1=a3m.readlines()
    a3m1=[x.strip() for x in a3m1]
with open(msa2, 'r') as a3m:
    a3m2=a3m.readlines()
    a3m2=[x.strip() for x in a3m2]


new_a3m1=[a3m1[0],a3m1[1]]
new_a3m2=[a3m2[0],a3m2[1]]

for i in range(len(order1)):
    #make sure that jackhmmer a3m header exist in the homolog file.
    name_a3m1 = [m for m in a3m1 if order1[i] in m]
    name_a3m2 = [n for n in a3m2 if order2[i] in n]
    ##make sure the ith sequence appears in both msa.
    if name_a3m1 and name_a3m2:
        name_pos1 = a3m1.index(name_a3m1[0])  #
        name_pos2 = a3m2.index(name_a3m2[0])
        seq1 = a3m1[name_pos1 + 1]
        seq2 = a3m2[name_pos2 + 1]

        new_a3m1.extend([name_a3m1[0], seq1])
        new_a3m2.extend([name_a3m2[0], seq2])


    elif name_a3m1 == [] and name_a3m2 != []:
        # remove the corresponding line in name_a3m2
        # get the position and delete the name and the seq from a3m file
        del_pos = a3m2.index(name_a3m2[0])
        del a3m2[del_pos]
        del a3m2[del_pos]
        print('case1:' + str(i) + ', remove '+ name_a3m2[0] +' from '+msa2)
        # should remove the nth seq from both fasta file and msa.
        # fasta
        del1_posA = fasta1.index([m for m in fasta1 if order1[i] in m][0])
        del1_posB = del1_posA + 2        
 
        del2_posA = fasta2.index([m for m in fasta2 if order2[i] in m][0])
        del2_posB = del2_posA + 2

        del fasta1[del1_posA:del1_posB]
        del fasta2[del2_posA:del2_posB]

    elif name_a3m2 == [] and name_a3m1 != []:

        del_pos = a3m1.index(name_a3m1[0])
        del a3m1[del_pos]
        del a3m1[del_pos]
        print("case2:" + str(i) + ', remove '+ name_a3m1[0] +' from '+msa1)
        del1_posA = fasta1.index([m for m in fasta1 if order1[i] in m][0])
        del1_posB = del1_posA + 2

        del2_posA = fasta2.index([m for m in fasta2 if order2[i] in m][0])
        del2_posB = del2_posA + 2

        del fasta1[del1_posA:del1_posB]
        del fasta2[del2_posA:del2_posB]
    else:
        print("case3:" + str(i) + ', remove both seqs '+order1[i]+' and '+order2[i])
        pass

#print(len(new_a3m1),len(new_a3m1))


#overwrite to the new file
with open(pr1, 'w+') as f:
    for item in fasta1:
        f.write("%s\n" % item)
with open(pr2, 'w+') as f:
    for item in fasta2:
        f.write("%s\n" % item)
with open(msa1, 'w+') as f:
    for item in new_a3m1:
        f.write("%s\n" % item)
with open(msa2, 'w+') as f:
    for item in new_a3m2:
        f.write("%s\n" % item)



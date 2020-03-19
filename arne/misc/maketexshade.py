#!/usr/bin/env python3

TAIL="""
\end{document}
"""
HEAD="""
\\documentclass[12pt]{article}

\\usepackage{texshade}


\\definecolor{colA}{rgb}{0.63,0,0.29} % A1004A
\\definecolor{colB}{rgb}{1.,0.776,0.} % FFC600
\\definecolor{colC}{rgb}{0.843,0.659,0.} % D7A800
\\definecolor{colD}{rgb}{0.592,0.447,0.} % 977200
\\definecolor{colE}{rgb}{0.357,0.408,0.} % 5B6800
\\definecolor{colF}{rgb}{0.357,0.408,0.} % 5B6800
\\definecolor{colG}{rgb}{0.357,0.408,0.} % 5B6800
\\definecolor{colH}{rgb}{0.373,0.541,0.} % 5E8A00
\\definecolor{colI}{rgb}{0.302,0.749,0.} % 43BF00
\\definecolor{colJ}{rgb}{0.,1.,0.} % 00FF00

%\\newcommand{\\colA}{ultramarine}
%\\newcommand{\\colB}{RubineRed}
%\\newcommand{\\colC}{Salmon}
%\\newcommand{\\colD}{VioletRed}
%\\newcommand{\\colE}{Mulberry}
%\\newcommand{\\colF}{Fuchsia}
%\\newcommand{\\colG}{DarkOrchid}
%\\newcommand{\\colH}{RoyalPurple}
%\\newcommand{\\colI}{CornflowerBlue}
%\\newcommand{\\colJ}{Blue}



\\newcommand{\\colA}{colA}
\\newcommand{\\colB}{colB}
\\newcommand{\\colC}{colC}
\\newcommand{\\colD}{colD}
\\newcommand{\\colE}{colE}
\\newcommand{\\colF}{colF}
\\newcommand{\\colG}{colG}
\\newcommand{\\colH}{colH}
\\newcommand{\\colI}{colI}
\\newcommand{\\colJ}{colJ}



\\begin{document}


"""

import os
import re
import argparse
from argparse import RawTextHelpFormatter

import os.path
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description =
                                 '- Makes an latexfile to be used with texshade for coloring an msa -',
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-color', required= True, help='Fasta file with sequence codes as numbers')
parser.add_argument('-MSA', required= True, help='MSA file in any format that texshade can read')
ns = parser.parse_args()
colorfile = ns.color
msafile = ns.MSA

# Open MSA filein fas
handle = open(colorfile, 'rU')


mapping={"1":"\colA","2":"\colB","3":"\colC","4":"\colD","5":"\colE","6":"\colF","7":"\colG","8":"\colH","9":"\colI"}


print (HEAD)
print ("\\begin{texshade}{"+msafile+"}")

seqnum=0
for record in SeqIO.parse(handle, 'fasta') :
    #print (record.name)
    seqnum+=1
    j=0
    for i in record.seq:
        j+=1
        print ("\\shaderegion{"+record.name+"}{"+str(j)+".."+str(j)+"}{"+mapping[i]+"}{White}")



print ("\\hideconsensus")
print ("\\end{texshade}")

print (TAIL)

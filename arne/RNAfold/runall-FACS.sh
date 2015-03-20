#!/bin/bash -x

# Extract Sequences form genebank file (hardcoded in script)
python bin/extract-nucleotide.py

# Experimemntal dat found in CSV files from Dan Daley

# Generate all mfold data etc
bin/FACS-analysis.pl araH    > araH-FACS.txt 
bin/FACS-analysis.pl narK    > narK-FACS.txt 

# Make plots

data <- read.csv(file="araH-FACS.txt",header=FALSE,sep=" ");

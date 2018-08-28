#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from Bio import SwissProt
from Bio import SwissProt
for record in SwissProt.parse(open("uniprot_sprot.dat")):
        print record.gene_name

from Bio import ExPASy
handle = ExPASy.get_sprot_raw(myaccessionnumber)
        


Holds information from a SwissProt record.

Members:
entry_name        Name of this entry, e.g. RL1_ECOLI.
data_class        Either 'STANDARD' or 'PRELIMINARY'.
molecule_type     Type of molecule, 'PRT',
sequence_length   Number of residues.

accessions        List of the accession numbers, e.g. ['P00321']
created           A tuple of (date, release).
sequence_update   A tuple of (date, release).
annotation_update A tuple of (date, release).

description       Free-format description.
gene_name         Gene name.  See userman.txt for description.
organism          The source of the sequence.
organelle         The origin of the sequence.
organism_classification  The taxonomy classification.  List of strings.
                         (http://www.ncbi.nlm.nih.gov/Taxonomy/)
                         taxonomy_id       A list of NCBI taxonomy id's.
host_organism     A list of names of the hosts of a virus, if any.
host_taxonomy_id  A list of NCBI taxonomy id's of the hosts, if any.
                         references        List of Reference objects.
                         comments          List of strings.
                         cross_references  List of tuples (db, id1[, id2][, id3]).  See the docs.
                         keywords          List of the keywords.
                         features          List of tuples (key name, from, to, description).
                                           from and to can be either integers for the residue
                                                             numbers, '<', '>', or '?'

                                                             seqinfo           tuple of (length, molecular weight, CRC32 value)
                                                             sequence          The sequence.
                                                             

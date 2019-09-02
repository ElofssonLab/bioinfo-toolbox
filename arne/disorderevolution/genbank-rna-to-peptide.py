#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import re
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC
import subprocess
from Bio import SwissProt

record = next(SeqIO.parse(filename, "genbank"))
from Bio import SeqIO
for seq_record in SeqIO.parse(filename, "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

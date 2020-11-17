#!/usr/bin/env python3
"""Docstring"""
import sys
#sys.path.append("/apps/pyconfold/")
import pyconsFold

fa1_file = sys.argv[1]
fa2_file = sys.argv[2]
rr_file = sys.argv[3]
out_folder = sys.argv[4]

pyconsFold.model_dock(fa1_file,fa2_file, out_folder, rr=rr_file)


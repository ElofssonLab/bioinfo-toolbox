#!/usr/bin/env python3
"""Docstring"""
import sys
sys.path.append("/apps/pyconfold/")
import pyconfold

fa_file = sys.argv[1]
rr_file = sys.argv[2]
out_folder = sys.argv[3]

pyconfold.model_dist(fa_file, rr_file, out_folder)


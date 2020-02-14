#!/usr/bin/env python3
import h5py
import math
import copy
import random
import datetime as day
import numpy as np
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
from Bio import AlignIO
from numpy import genfromtxt
from keras import backend as K
from keras.regularizers import l2
from keras.layers.merge import add, subtract, average, multiply
from keras.layers import *

ATOMIC_WEIGHTS = {'H':1.008, 'HE':4.002602, 'LI':6.94, 'BE':9.012182,
       'B':10.81, 'C':12.011, 'N':14.007, 'O':15.999, 'F':18.9984032,
       'NE':20.1797, 'NA':22.98976928, 'MG':24.305, 'AL':26.9815386,
       'SI':28.085, 'P':30.973762, 'S':32.06, 'CL':35.45, 'AR':39.948,
       'K':39.0983, 'CA':40.078, 'SC':44.955912, 'TI':47.867, 'V':50.9415,
       'CR':51.9961, 'MN':54.938045, 'FE':55.845, 'CO':58.933195,
       'NI':58.6934, 'CU':63.546, 'ZN':65.38, 'GA':69.723, 'GE':72.630,
       'AS':74.92160, 'SE':78.96, 'BR':79.904, 'RB':85.4678, 'SR':87.62,
       'Y':88.90585, 'ZR':91.224, 'NB':92.90638, 'MO':95.96, 'TC':98,
       'RU':101.07, 'RH':102.90550, 'PD':106.42, 'AG':107.8682, 'CD':112.411,
       'IN':114.818, 'SN':118.710, 'SB':121.760, 'TE':127.60, 'I':126.90447,
       'XE':131.293, 'CS':132.9054519, 'BA':137.327, 'LA':138.90547,
       'CE':140.116, 'PR':140.90765, 'ND':144.242, 'PM':145, 'SM':150.36,
       'EU':151.964, 'GD':157.25, 'TB':158.92535, 'DY':162.500, 'HO':164.93032,
       'ER':167.259, 'TM':168.93421, 'YB':173.054, 'LU':174.9668, 'HF':178.49,
       'TA':180.94788, 'W':183.84, 'RE':186.207, 'OS':190.23, 'IR':192.217,
       'PT':195.084, 'AU':196.966569, 'HG':200.592, 'TL':204.38, 'PB':207.2,
       'BI':208.98040, 'PO':209, 'AT':210, 'RN':222, 'FR':223, 'RA':226,
       'AC':227, 'TH':232.03806, 'PA':231.03588, 'U':238.02891, 'NP':237,
       'PU':244, 'AM':243, 'CM':247, 'BK':247, 'CF':251, 'ES':252, 'FM':257,
       'MD':258, 'NO':259, 'LR':262, 'RF':267, 'DB':268, 'SG':269, 'BH':270,
       'HS':269, 'MT':278, 'DS':281, 'RG':281, 'CN':285, 'UUT':286, 'FL':289,
       'UUP':288, 'LV':293, 'UUS':294}

'''
ID Type Description
1 Sulphur/selenium CYS:SG, MET:SD, MSE:SE
2 Nitrogen (amide) ASN:ND2, GLN:NE2, backbone N (including N-terminal)
3 Nitrogen (aromatic) HIS:ND1/NE2, TRP:NE1
4 Nitrogen (guanidinium) ARG:NE/NH*
5 Nitrogen (ammonium) LYS:NZ
6 Oxygen (carbonyl) ASN:OD1, GLN:OE1, backbone O (except C-terminal)
7 Oxygen (hydroxyl) SER:OG, THR:OG1, TYR:OH
8 Oxygen (carboxyl) ASP:OD*, GLU:OE*, C-terminal O/OXT
9 Carbon (sp2) ARG:CZ, ASN:CG, ASP:CG, GLN:CD, GLU:CD, backbone C
10 Carbon (aromatic) HIS:CG/CD2/CE1, PHE:CG/CD*/CE*/CZ,
TRP:CG/CD*/CE*/CZ*/CH2, TYR:CG/CD*/CE*/CZ
11 Carbon (sp3) ALA:CB, ARG:CB/CG/CD, ASN:CB, ASP:CB, CYS:CB,
GLN:CB/CG, GLU:CB/CG, HIS:CB, ILE:CB/CG*/CD1,
LEU:CB/CG/CD*, LYS:CB/CG/CD/CE, MET:CB/CG/CE,
MSE:CB/CG/CE, PHE:CB, PRO:CB/CG/CD, SER:CB,
THR:CB/CG2, TRP:CB, TYR:CB, VAL:CB/CG*,
'''

allowed_atoms = ['N','CA','C','O','OXT','CB',
                 'CZ','CZ1','CZ2','CE','CE1','CE2',
                 'CD','CD1','CD2','CG','CG1','CG2','CG3',
                 'CH2','OE','OE1','OE2','OD','OD1','OD2',
                 'OH','OG','OG1','OG2','NZ','NE','NE1','NE2',
                 'NH','NH1','NH2','ND','ND1','ND2',
                 'SG','SD','SE']

atom_types = {'ALA':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'ARG':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'NE':[0,0,0,1,0,0,0,0,0,0,0],
                     'NH1':[0,0,0,1,0,0,0,0,0,0,0],
                     'NH2':[0,0,0,1,0,0,0,0,0,0,0],
                      'CZ':[0,0,0,0,0,0,0,0,1,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CD':[0,0,0,0,0,0,0,0,0,0,1]},
              'ASN':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                     'ND2':[0,1,0,0,0,0,0,0,0,0,0], 
                     'OD1':[0,0,0,0,0,1,0,0,0,0,0], 
                      'CG':[0,0,0,0,0,0,0,0,1,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'ASP':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                     'OD1':[0,0,0,0,0,0,0,1,0,0,0], 
                     'OD2':[0,0,0,0,0,0,0,1,0,0,0],
                      'CG':[0,0,0,0,0,0,0,0,1,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'CYS':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'SG':[1,0,0,0,0,0,0,0,0,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'GLN':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                     'NE2':[0,1,0,0,0,0,0,0,0,0,0], 
                     'OE1':[0,0,0,0,0,1,0,0,0,0,0], 
                      'CD':[0,0,0,0,0,0,0,0,1,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,0,1]},
              'GLU':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                     'OE1':[0,0,0,0,0,0,0,1,0,0,0],
                     'OE2':[0,0,0,0,0,0,0,1,0,0,0],
                      'CD':[0,0,0,0,0,0,0,0,1,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1],
                      'CG':[0,0,0,0,0,0,0,0,0,0,1]},
              'GLY':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0]},
              'HIS':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                     'ND1':[0,0,1,0,0,0,0,0,0,0,0], 
                     'NE2':[0,0,1,0,0,0,0,0,0,0,0], 
                      'CG':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CD2':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CE1':[0,0,0,0,0,0,0,0,0,1,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'ILE':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                     'CG1':[0,0,0,0,0,0,0,0,0,0,1], 
                     'CG2':[0,0,0,0,0,0,0,0,0,0,1],
                     'CD1':[0,0,0,0,0,0,0,0,0,0,1]},
              'LEU':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,0,1], 
                     'CD1':[0,0,0,0,0,0,0,0,0,0,1],
                     'CD2':[0,0,0,0,0,0,0,0,0,0,1]},
              'LYS':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'NZ':[0,0,0,0,1,0,0,0,0,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CD':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CE':[0,0,0,0,0,0,0,0,0,0,1]},
              'MET':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CE':[0,0,0,0,0,0,0,0,0,0,1], 
                      'SD':[1,0,0,0,0,0,0,0,0,0,0]},
              'PHE':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CD1':[0,0,0,0,0,0,0,0,0,1,0],
                     'CD2':[0,0,0,0,0,0,0,0,0,1,0],
                     'CE1':[0,0,0,0,0,0,0,0,0,1,0],
                     'CE2':[0,0,0,0,0,0,0,0,0,1,0],
                      'CZ':[0,0,0,0,0,0,0,0,0,1,0]},
              'PRO':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CD':[0,0,0,0,0,0,0,0,0,0,1]},
              'SER':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'OG':[0,0,0,0,0,0,1,0,0,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'THR':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                     'OG1':[0,0,0,0,0,0,1,0,0,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1],
                     'CG2':[0,0,0,0,0,0,0,0,0,0,1]},
              'TRP':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                      'CG':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CD1':[0,0,0,0,0,0,0,0,0,1,0],
                     'CD2':[0,0,0,0,0,0,0,0,0,1,0],
                     'CE1':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CE2':[0,0,0,0,0,0,0,0,0,1,0],
                     'CE3':[0,0,0,0,0,0,0,0,0,1,0],
                     'CZ2':[0,0,0,0,0,0,0,0,0,1,0],
                     'CZ3':[0,0,0,0,0,0,0,0,0,1,0],
                     'CH2':[0,0,0,0,0,0,0,0,0,1,0],
                     'NE1':[0,0,1,0,0,0,0,0,0,0,0]},
              'TYR':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CG':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CD1':[0,0,0,0,0,0,0,0,0,1,0],
                     'CD2':[0,0,0,0,0,0,0,0,0,1,0],
                     'CE1':[0,0,0,0,0,0,0,0,0,1,0], 
                     'CE2':[0,0,0,0,0,0,0,0,0,1,0],
                      'CZ':[0,0,0,0,0,0,0,0,0,1,0], 
                      'OH':[0,0,0,0,0,0,1,0,0,0,0], 
                      'CB':[0,0,0,0,0,0,0,0,0,0,1]},
              'VAL':{  'N':[0,1,0,0,0,0,0,0,0,0,0],
                      'CA':[0,0,0,0,0,0,0,0,1,0,0],
                       'C':[0,0,0,0,0,0,0,0,1,0,0],
                       'O':[0,0,0,0,0,1,0,0,0,0,0],
                     'OXT':[0,0,0,0,0,0,0,1,0,0,0],
                      'CB':[0,0,0,0,0,0,0,0,0,0,1], 
                     'CG1':[0,0,0,0,0,0,0,0,0,0,1],
                     'CG2':[0,0,0,0,0,0,0,0,0,0,1]}}


max_asa = {'A':121.0,
           'R':265.0,
           'N':187.0,
           'D':187.0,
           'C':148.0,
           'Q':214.0,
           'E':214.0,
           'G':97.0,
           'H':216.0,
           'I':195.0,
           'L':191.0,
           'K':230.0,
           'M':203.0,
           'F':228.0,
           'P':154.0,
           'S':143.0,
           'T':163.0,
           'W':264.0,
           'Y':255.0,
           'V':165.0,
           'X':265.0,
           'U':265.0,
           'O':265.0}



atom_type = {'C':[1.0,0.0,0.0,0.0],
             'O':[0.0,1.0,0.0,0.0],
             'N':[0.0,0.0,1.0,0.0],
             'S':[0.0,0.0,0.0,1.0]}

onechar = ['A','R','N','D','C','Q','E','G','H','I',
           'L','K','M','F','P','S','T','W','Y','V']

threechar = ['ALA','ARG','ASN','ASP',
             'CYS','GLN','GLU','GLY',
             'HIS','ILE','LEU','LYS',
             'MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']

one2three = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP',
             'C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
             'H':'HIS','I':'ILE','L':'LEU','K':'LYS',
             'M':'MET','F':'PHE','P':'PRO','S':'SER',
             'T':'THR','W':'TRP','Y':'TYR','V':'VAL',
             'X':'XXX','U':'XXX','O':'XXX'}

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'XXX':'X'}

aa_encode = {
'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'R':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'N':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'D':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'C':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'Q':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'E':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
'G':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
'H':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
'I':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
'L':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
'K':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
'M':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
'F':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
'P':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
'X':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'U':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'O':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}

embed7d = {
'A':[1.00, 0.19, 0.25, 0.56, 1.00, 0.86, 0.19],
'R':[0.00, 0.76, 0.66, 0.34, 0.00, 1.00, 0.01],
'N':[0.17, 0.27, 0.68, 0.19, 0.18, 0.68, 0.18],
'D':[0.54, 0.29, 0.22, 0.51, 0.41, 0.11, 0.00],
'C':[0.88, 0.03, 0.65, 0.09, 0.54, 0.65, 0.74],
'Q':[0.71, 0.69, 0.65, 0.46, 0.49, 0.84, 0.21],
'E':[0.68, 0.57, 0.37, 0.90, 0.87, 0.32, 0.01],
'G':[0.90, 0.00, 0.00, 0.05, 0.70, 0.59, 0.17],
'H':[0.20, 0.34, 0.92, 0.39, 0.25, 0.76, 0.39],
'I':[0.71, 0.45, 0.89, 0.94, 0.83, 0.35, 1.00],
'L':[0.90, 0.76, 0.79, 0.70, 0.78, 0.89, 0.49],
'K':[0.15, 0.42, 0.86, 0.27, 0.11, 0.86, 0.08],
'M':[0.88, 0.78, 0.79, 0.66, 0.67, 0.89, 0.42],
'F':[0.83, 0.77, 0.88, 0.69, 0.50, 0.81, 0.58],
'P':[0.24, 0.05, 0.82, 0.00, 0.82, 0.00, 0.99],
'S':[0.76, 0.11, 0.37, 0.19, 0.55, 0.57, 0.41],
'T':[0.37, 0.18, 0.73, 0.23, 0.55, 0.38, 0.67],
'W':[0.95, 1.00, 1.00, 1.00, 0.37, 0.89, 0.66],
'Y':[0.66, 0.67, 0.89, 0.58, 0.23, 0.70, 0.52],
'V':[0.71, 0.18, 0.79, 0.48, 0.77, 0.30, 0.99],
'X':[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
'U':[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
'O':[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]}

ss3_encode = {'H':[1.0,0.0,0.0],
              'G':[1.0,0.0,0.0],
              'I':[1.0,0.0,0.0],
              'E':[0.0,1.0,0.0],
              'B':[0.0,1.0,0.0],
              'T':[0.0,0.0,1.0],
              'S':[0.0,0.0,1.0],
              'C':[0.0,0.0,1.0]}

ss8_encode = {'H':[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
              'G':[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],
              'I':[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0],
              'E':[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],
              'B':[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],
              'T':[0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0],
              'S':[0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],
              'C':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]}

mdb_encode = {'S':0.0, 'D':1.0, 'C':0.0, 'X':0.5}

#thrlist = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0]

thrlist = list(np.arange(0.0,1.0,0.01))

blosum62 = {'ALA':{'ALA':0.2901,'ARG':0.031,'ASN':0.0256,'ASP':0.0297,'CYS':0.0216,'GLN':0.0256,'GLU':0.0405,'GLY':0.0783,'HIS':0.0148,'ILE':0.0432,'LEU':0.0594,'LYS':0.0445,'MET':0.0175,'PHE':0.0216,'PRO':0.0297,'SER':0.085,'THR':0.0499,'TRP':0.0054,'TYR':0.0175,'VAL':0.0688},
'ARG':{'ALA':0.0446,'ARG':0.345,'ASN':0.0388,'ASP':0.031,'CYS':0.0078,'GLN':0.0484,'GLU':0.0523,'GLY':0.0329,'HIS':0.0233,'ILE':0.0233,'LEU':0.0465,'LYS':0.1202,'MET':0.0155,'PHE':0.0174,'PRO':0.0194,'SER':0.0446,'THR':0.0349,'TRP':0.0058,'TYR':0.0174,'VAL':0.031},
'ASN':{'ALA':0.0427,'ARG':0.0449,'ASN':0.3169,'ASP':0.0831,'CYS':0.009,'GLN':0.0337,'GLU':0.0494,'GLY':0.0652,'HIS':0.0315,'ILE':0.0225,'LEU':0.0315,'LYS':0.0539,'MET':0.0112,'PHE':0.018,'PRO':0.0202,'SER':0.0697,'THR':0.0494,'TRP':0.0045,'TYR':0.0157,'VAL':0.027},
'ASP':{'ALA':0.041,'ARG':0.0299,'ASN':0.069,'ASP':0.3974,'CYS':0.0075,'GLN':0.0299,'GLU':0.0914,'GLY':0.0466,'HIS':0.0187,'ILE':0.0224,'LEU':0.028,'LYS':0.0448,'MET':0.0093,'PHE':0.0149,'PRO':0.0224,'SER':0.0522,'THR':0.0354,'TRP':0.0037,'TYR':0.0112,'VAL':0.0243},
'CYS':{'ALA':0.065,'ARG':0.0163,'ASN':0.0163,'ASP':0.0163,'CYS':0.4837,'GLN':0.0122,'GLU':0.0163,'GLY':0.0325,'HIS':0.0081,'ILE':0.0447,'LEU':0.065,'LYS':0.0203,'MET':0.0163,'PHE':0.0203,'PRO':0.0163,'SER':0.0407,'THR':0.0366,'TRP':0.0041,'TYR':0.0122,'VAL':0.0569},
'GLN':{'ALA':0.0559,'ARG':0.0735,'ASN':0.0441,'ASP':0.0471,'CYS':0.0088,'GLN':0.2147,'GLU':0.1029,'GLY':0.0412,'HIS':0.0294,'ILE':0.0265,'LEU':0.0471,'LYS':0.0912,'MET':0.0206,'PHE':0.0147,'PRO':0.0235,'SER':0.0559,'THR':0.0412,'TRP':0.0059,'TYR':0.0206,'VAL':0.0353},
'GLU':{'ALA':0.0552,'ARG':0.0497,'ASN':0.0405,'ASP':0.0902,'CYS':0.0074,'GLN':0.0645,'GLU':0.2965,'GLY':0.035,'HIS':0.0258,'ILE':0.0221,'LEU':0.0368,'LYS':0.0755,'MET':0.0129,'PHE':0.0166,'PRO':0.0258,'SER':0.0552,'THR':0.0368,'TRP':0.0055,'TYR':0.0166,'VAL':0.0313},
'GLY':{'ALA':0.0783,'ARG':0.0229,'ASN':0.0391,'ASP':0.0337,'CYS':0.0108,'GLN':0.0189,'GLU':0.0256,'GLY':0.5101,'HIS':0.0135,'ILE':0.0189,'LEU':0.0283,'LYS':0.0337,'MET':0.0094,'PHE':0.0162,'PRO':0.0189,'SER':0.0513,'THR':0.0297,'TRP':0.0054,'TYR':0.0108,'VAL':0.0243},
'HIS':{'ALA':0.042,'ARG':0.0458,'ASN':0.0534,'ASP':0.0382,'CYS':0.0076,'GLN':0.0382,'GLU':0.0534,'GLY':0.0382,'HIS':0.355,'ILE':0.0229,'LEU':0.0382,'LYS':0.0458,'MET':0.0153,'PHE':0.0305,'PRO':0.0191,'SER':0.042,'THR':0.0267,'TRP':0.0076,'TYR':0.0573,'VAL':0.0229},
'ILE':{'ALA':0.0471,'ARG':0.0177,'ASN':0.0147,'ASP':0.0177,'CYS':0.0162,'GLN':0.0133,'GLU':0.0177,'GLY':0.0206,'HIS':0.0088,'ILE':0.271,'LEU':0.1679,'LYS':0.0236,'MET':0.0368,'PHE':0.0442,'PRO':0.0147,'SER':0.025,'THR':0.0398,'TRP':0.0059,'TYR':0.0206,'VAL':0.1767},
'LEU':{'ALA':0.0445,'ARG':0.0243,'ASN':0.0142,'ASP':0.0152,'CYS':0.0162,'GLN':0.0162,'GLU':0.0202,'GLY':0.0213,'HIS':0.0101,'ILE':0.1154,'LEU':0.3755,'LYS':0.0253,'MET':0.0496,'PHE':0.0547,'PRO':0.0142,'SER':0.0243,'THR':0.0334,'TRP':0.0071,'TYR':0.0223,'VAL':0.0962},
'LYS':{'ALA':0.057,'ARG':0.1071,'ASN':0.0415,'ASP':0.0415,'CYS':0.0086,'GLN':0.0535,'GLU':0.0708,'GLY':0.0432,'HIS':0.0207,'ILE':0.0276,'LEU':0.0432,'LYS':0.2781,'MET':0.0155,'PHE':0.0155,'PRO':0.0276,'SER':0.0535,'THR':0.0397,'TRP':0.0052,'TYR':0.0173,'VAL':0.0328},
'MET':{'ALA':0.0522,'ARG':0.0321,'ASN':0.0201,'ASP':0.0201,'CYS':0.0161,'GLN':0.0281,'GLU':0.0281,'GLY':0.0281,'HIS':0.0161,'ILE':0.1004,'LEU':0.1968,'LYS':0.0361,'MET':0.1606,'PHE':0.0482,'PRO':0.0161,'SER':0.0361,'THR':0.0402,'TRP':0.008,'TYR':0.0241,'VAL':0.0924},
'PHE':{'ALA':0.0338,'ARG':0.019,'ASN':0.0169,'ASP':0.0169,'CYS':0.0106,'GLN':0.0106,'GLU':0.019,'GLY':0.0254,'HIS':0.0169,'ILE':0.0634,'LEU':0.1142,'LYS':0.019,'MET':0.0254,'PHE':0.3869,'PRO':0.0106,'SER':0.0254,'THR':0.0254,'TRP':0.0169,'TYR':0.0888,'VAL':0.055},
'PRO':{'ALA':0.0568,'ARG':0.0258,'ASN':0.0233,'ASP':0.031,'CYS':0.0103,'GLN':0.0207,'GLU':0.0362,'GLY':0.0362,'HIS':0.0129,'ILE':0.0258,'LEU':0.0362,'LYS':0.0413,'MET':0.0103,'PHE':0.0129,'PRO':0.4935,'SER':0.0439,'THR':0.0362,'TRP':0.0026,'TYR':0.0129,'VAL':0.031},
'SER':{'ALA':0.1099,'ARG':0.0401,'ASN':0.0541,'ASP':0.0489,'CYS':0.0175,'GLN':0.0332,'GLU':0.0524,'GLY':0.0663,'HIS':0.0192,'ILE':0.0297,'LEU':0.0419,'LYS':0.0541,'MET':0.0157,'PHE':0.0209,'PRO':0.0297,'SER':0.2199,'THR':0.082,'TRP':0.0052,'TYR':0.0175,'VAL':0.0419},
'THR':{'ALA':0.073,'ARG':0.0355,'ASN':0.0434,'ASP':0.0375,'CYS':0.0178,'GLN':0.0276,'GLU':0.0394,'GLY':0.0434,'HIS':0.0138,'ILE':0.0533,'LEU':0.0651,'LYS':0.0454,'MET':0.0197,'PHE':0.0237,'PRO':0.0276,'SER':0.0927,'THR':0.2465,'TRP':0.0059,'TYR':0.0178,'VAL':0.071},
'TRP':{'ALA':0.0303,'ARG':0.0227,'ASN':0.0152,'ASP':0.0152,'CYS':0.0076,'GLN':0.0152,'GLU':0.0227,'GLY':0.0303,'HIS':0.0152,'ILE':0.0303,'LEU':0.053,'LYS':0.0227,'MET':0.0152,'PHE':0.0606,'PRO':0.0076,'SER':0.0227,'THR':0.0227,'TRP':0.4924,'TYR':0.0682,'VAL':0.0303},
'TYR':{'ALA':0.0405,'ARG':0.028,'ASN':0.0218,'ASP':0.0187,'CYS':0.0093,'GLN':0.0218,'GLU':0.028,'GLY':0.0249,'HIS':0.0467,'ILE':0.0436,'LEU':0.0685,'LYS':0.0312,'MET':0.0187,'PHE':0.1308,'PRO':0.0156,'SER':0.0312,'THR':0.028,'TRP':0.028,'TYR':0.3178,'VAL':0.0467},
'VAL':{'ALA':0.07,'ARG':0.0219,'ASN':0.0165,'ASP':0.0178,'CYS':0.0192,'GLN':0.0165,'GLU':0.0233,'GLY':0.0247,'HIS':0.0082,'ILE':0.1646,'LEU':0.1303,'LYS':0.0261,'MET':0.0316,'PHE':0.0357,'PRO':0.0165,'SER':0.0329,'THR':0.0494,'TRP':0.0055,'TYR':0.0206,'VAL':0.2689}}

def file_load(filename):
##### given a path, opens the corresponding file or print out it is 'not found!'

    try: infile = open(filename,'r')
    except:
        print (filename+' not found!')
        return 'none'
    return infile

def get_fasta(code, chain):

    line = ''
    seq = '>'+code+'\n'
    for residue in chain:
        if not is_aa(residue): continue
        line += three2one[residue.get_resname()]
        if len(line) == 80:
            line += '\n'
            seq += line
            line = ''
    if line != '': 
        seq += line
    return seq

def get_interactions(chain1, chain2, thr):

    ppi_residues = {'chain1':[], 'chain2':[]}
    dssp = DSSP(model, "/local-pdb/1mot.pdb")

    atoms_chain1 = []
    for residue1 in chain1:
        if not is_aa(residue1): continue
        for atom1 in residue1: atoms_chain1.append(atom1)

    for residue2 in chain2:
        if not is_aa(residue2): continue
        for atom2 in residue2: residue2_contacts = NeighborSearch(atoms_chain1).search(atom2.get_coord(), thr, level='R')
        if residue2_contacts != []: ppi_residues['chain2'].append(residue2)
        for residue1 in residue2_contacts: 
            if residue1 not in ppi_residues['chain1']: ppi_residues['chain1'].append(residue1)
            
    return ppi_residues

def center_of_mass(chain):
    masses = []
    coordinates = []
    center = [None, None, None]

    for residue in chain:
        for atom in residue:
            coordinates.append([atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]])
            try:
                element_name = atom.get_name()[0].lstrip('1234567890')
                masses.append(ATOMIC_WEIGHTS[element_name])
            except: 
                element_name = atom.get_name()[0:2].lstrip('1234567890')
                masses.append(ATOMIC_WEIGHTS[element_name])

    total_mass = sum(masses)
    weights = [float(atom_mass)/total_mass for atom_mass in masses]

    center = [sum([coordinates[i][j] * weights[i] 
          for i in range(len(weights))]) for j in range(3)]
    center_rounded = [round(center[i], 3) for i in range(3)]

    return center_rounded

def interaction_center_of_mass(ppi_residues):
    masses = []
    coordinates = []
    center = [None, None, None]

    for residue in ppi_residues['chain1']:
        for atom in residue:
            coordinates.append([atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]])
            try:
                element_name = atom.get_name()[0].lstrip('1234567890')
                masses.append(ATOMIC_WEIGHTS[element_name])
            except:
                element_name = atom.get_name()[0:2].lstrip('1234567890')
                masses.append(ATOMIC_WEIGHTS[element_name])

    for residue in ppi_residues['chain2']:
        for atom in residue:
            coordinates.append([atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]])
            try:
                element_name = atom.get_name()[0].lstrip('1234567890')
                masses.append(ATOMIC_WEIGHTS[element_name])
            except:
                element_name = atom.get_name()[0:2].lstrip('1234567890')
                masses.append(ATOMIC_WEIGHTS[element_name])

    total_mass = sum(masses)
    weights = [float(atom_mass)/total_mass for atom_mass in masses]

    center = [sum([coordinates[i][j] * weights[i]
          for i in range(len(weights))]) for j in range(3)]
    center_rounded = [round(center[i], 3) for i in range(3)]
    return center_rounded

def midpoint_find(seta, setb):
##### returns the coordinates of the midpoint between two points, given their 3D coordinates #####
    midpoint = []
    for ax in range(len(seta)): midpoint.append((seta[ax]+setb[ax])/2)
    return midpoint


def distance_find(seta, setb):
##### returns the distance between two points, given their 3D coordinates #####
    dab = math.sqrt((seta[0]-setb[0])**2+(seta[1]-setb[1])**2+(seta[2]-setb[2])**2)
    return dab

def rotation_find(c1, c2):
##### returns the rotation around c1 that brings c2 on the positive Y axis #####
    v = [c2[0]-c1[0], c2[1]-c1[1], c2[2]-c1[2]]
    x = [1.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [0.0, 0.0, 1.0]

    dotvz = v[0]*z[0] + v[1]*z[1] + v[2]*z[2]
    detvz = v[0]*z[1]*y[2] + z[0]*y[1]*v[2] + y[0]*v[1]*z[2] - v[2]*z[1]*y[0] - z[2]*y[1]*v[0] - y[2]*v[1]*z[0]
    angvz = math.atan2(detvz, dotvz)

    if angvz < 0: angvz = math.radians(360.0)+angvz
    yrot = np.array(([[math.cos(angvz), 0, math.sin(angvz)],[0,1,0],[-(math.sin(angvz)),0,math.cos(angvz)]]), dtype=np.float64)

    v1 = yrot.dot(v)
    
    dotv1y = v1[0]*y[0] + v1[1]*y[1] + v1[2]*y[2] 
    detv1y = v1[0]*y[1]*x[2] + y[0]*x[1]*v1[2] + x[0]*v1[1]*y[2] - v1[2]*y[1]*x[0] - y[2]*x[1]*v1[0] - x[2]*v1[1]*y[0]
    angv1y = math.atan2(detv1y, dotv1y)

    if angv1y < 0: angv1y = math.radians(360.0)+angv1y
    xrot = np.array(([[1,0,0],[0,math.cos(angv1y),-(math.sin(angv1y))],[0,math.sin(angv1y),math.cos(angv1y)]]), dtype=np.float64)
    
    yxrot = np.dot(xrot, yrot)
    return yxrot

def interface_extract(chain1, chain2, thr, dssp1, dssp2):

    ppi_residues = get_interactions(chain1, chain2, thr)
    while len(ppi_residues['chain1']) == 0 or len(ppi_residues['chain2']) == 0:
        thr += 2
        ppi_residues = get_interactions(chain1, chain2, thr)

    #cm1 = np.array(center_of_mass(ppi_residues['chain1']))
    #cm2 = np.array(center_of_mass(ppi_residues['chain2']))
    #rot = rotation_find(cm1, cm2)

    #cm = np.array(interaction_center_of_mass(ppi_residues))
    #cm_r = np.dot(rot, (cm-cm1))

    ppi = {1:[], 2:[]}

    for residue1 in chain1: 
        ppi[1].append([residue1])
        if dssp1[(chain1.get_id(), residue1.get_id())][3] <= 0.3: continue 
        if residue1 in ppi_residues['chain1']: 
            residue_interaction = get_interactions([residue1], chain2, 12)
            for residue2 in residue_interaction['chain2']:
                caca = distance_find(residue1['CA'].get_coord(), residue2['CA'].get_coord())
                if residue1.get_resname() != 'GLY' and residue2.get_resname() != 'GLY': 
                    cbcb = distance_find(residue1['CB'].get_coord(), residue2['CB'].get_coord())
                elif residue1.get_resname() == 'GLY' and residue2.get_resname() == 'GLY':
                    cbcb = distance_find(residue1['CA'].get_coord(), residue2['CA'].get_coord())
                elif residue1.get_resname() == 'GLY': 
                    cbcb = distance_find(residue1['CA'].get_coord(), residue2['CB'].get_coord())
                elif residue2.get_resname() == 'GLY':
                    cbcb = distance_find(residue1['CB'].get_coord(), residue2['CA'].get_coord())
                ppi[1][-1].append([residue2, caca, cbcb])

    for residue2 in chain2:
        ppi[2].append([residue2])
        if dssp2[(chain2.get_id(), residue2.get_id())][3] <= 0.3: continue
        if residue2 in ppi_residues['chain2']:
            residue_interaction = get_interactions(chain1, [residue2], 12)
            for residue1 in residue_interaction['chain1']:
                caca = distance_find(residue1['CA'].get_coord(), residue2['CA'].get_coord())
                if residue2.get_resname() != 'GLY' and residue1.get_resname() != 'GLY': 
                    cbcb = distance_find(residue1['CB'].get_coord(), residue2['CB'].get_coord())
                elif residue2.get_resname() == 'GLY' and residue1.get_resname() == 'GLY':
                    cbcb = distance_find(residue2['CA'].get_coord(), residue1['CA'].get_coord())
                elif residue2.get_resname() == 'GLY':
                    cbcb = distance_find(residue2['CA'].get_coord(), residue1['CB'].get_coord())
                elif residue1.get_resname() == 'GLY':
                    cbcb = distance_find(residue2['CB'].get_coord(), residue1['CA'].get_coord())
                ppi[2][-1].append([residue1, caca, cbcb])
    return ppi

def box_extract(residue, own_chain, ppi, cube_side, cell_side=1, thr=3, features_number=13, influence=0):
    if own_chain == 1: opposite_chain = 2
    elif own_chain == 2: opposite_chain = 1

    for residue_group in ppi[own_chain]: 
        center_name = residue_group[0].get_resname()+str(residue_group[0].get_id()[1])
        if center_name == residue: 
            center_residue = residue_group
            break

    sub_interface = []
    while sub_interface == []:
        thr += 1
        for interacting_residue in center_residue[1:]:
            if interacting_residue[2] < thr: sub_interface.append(interacting_residue[0])
    sub_interface = []
    thr += 4
    for interacting_residue in center_residue[1:]:
        if interacting_residue[2] < thr: sub_interface.append(interacting_residue[0])

    residue_ca = np.array(center_residue[0]['CA'].get_coord())
    sub_interface_cm = np.array(center_of_mass(sub_interface))
    cube_center = np.array(midpoint_find(residue_ca, sub_interface_cm))
    rot = rotation_find(cube_center, sub_interface_cm)
    rotated_cm = np.dot(rot, np.array(sub_interface_cm-cube_center))

    half_diagonal = math.sqrt(cube_side**2+math.sqrt((cube_side**2)*2)**2)/2
    sphere = sphere_build(ppi[own_chain]+ppi[opposite_chain], cube_center, half_diagonal)

    cells = int(cube_side/cell_side)
    box = [[[[0.0 for f in range(features_number)] for x in range(cells)] for y in range(cells)] for z in range(cells)]
    box = np.array(box, dtype=np.float32)

    pdb = ''
    count = 0
    chain = 'A'
    for residue in sphere:
        name = residue.get_resname()
        for residue2 in ppi[opposite_chain]:
            if residue == residue2[0]:
                if chain == 'A': pdb += 'TER\n'
                chain = 'B'
                break

        if name not in atom_types: continue
        #position_scores = profile[residue_count]
        for atom in residue:
            if atom.get_name().lstrip('0123456789')[0]=='H': continue
            rotated_atom = np.dot(rot, (np.array(atom.get_coord()) - cube_center))

            xtop, xbot = range_find(rotated_atom[0], influence, cube_side)
            ytop, ybot = range_find(rotated_atom[1], influence, cube_side)
            ztop, zbot = range_find(rotated_atom[2], influence, cube_side)

            if xbot == -1 and xtop == -1: continue
            if ybot == -1 and ytop == -1: continue
            if zbot == -1 and ztop == -1: continue
            if xbot == 20 and xtop == 20: continue
            if ybot == 20 and ytop == 20: continue
            if zbot == 20 and ztop == 20: continue

            count += 1
            for x in range(xbot, xtop+1):
                for y in range(ybot, ytop+1):
                    for z in range(zbot, ztop+1):
                        cell_center = np.array([(x-9.5), (y-9.5), (z-9.5)])
                        #distance = distance_find(atom.get_coord(), cell_center)
                        for pos in range(len(atom_types[name][atom.get_name()])):
                            if atom_types[name][atom.get_name()][pos] == 1: box[x][y][z][pos] = 1.0 #+= math.exp(-((distance/dev)**2))

                        if chain == 'A': box[x][y][z][11] = 1.0
                        else: box[x][y][z][12] = 1.0

                        #if rotated_atom[0] < (x-10) or rotated_atom[0] > (x-9): continue
                        #if rotated_atom[1] < (y-10) or rotated_atom[1] > (y-9): continue
                        #if rotated_atom[2] < (z-10) or rotated_atom[2] > (z-9): continue

                        #pos = 13
                        #for code in threechar:
                        #    box[x][y][z][pos] = blosum62[name][code]

                        #pos = 13
                        #for code in onechar:
                        #    if code in position_scores: box[x][y][z][pos] = position_scores[code]
                        #    pos += 1
            #pdb += pdbline_format(count, atom.get_name(), rotated_atom, residue, chain)
    #pdb += 'TER\n'

    side_length = int(cube_side/cell_side)
    box = box.reshape(side_length, side_length, side_length, features_number)
    return box    

def sphere_build(chains, cube_center, half_diagonal):

    sphere = []
    for residue in chains:
        if residue[0].has_id('N'):
            if distance_find(cube_center, residue[0]['N'].get_coord()) <= half_diagonal: 
                sphere.append(residue[0])
                continue
        if residue[0].has_id('O'):
            if distance_find(cube_center, residue[0]['O'].get_coord()) <= half_diagonal: 
                sphere.append(residue[0])
                continue
        if residue[0].has_id('CB'):
            if distance_find(cube_center, residue[0]['CB'].get_coord()) <= half_diagonal: 
                sphere.append(residue[0])
                continue

    return sphere

def pdbline_format(atom_number, atom_name, coordinates, residue, chain):

    atom_number = str(atom_number)
    residue_name = residue.get_resname()
    residue_number = str(residue.get_id()[1])
    x = str(coordinates[0]).split('.')[0]+'.'+str(coordinates[0]).split('.')[1][:3]
    if 'e' in str(coordinates[0]):
        if int(str(coordinates[0]).split('e')[1]) < 0: x = '0.000'
    y = str(coordinates[1]).split('.')[0]+'.'+str(coordinates[1]).split('.')[1][:3]
    if 'e' in str(coordinates[1]):
        if int(str(coordinates[1]).split('e')[1]) < 0: y = '0.000'
    z = str(coordinates[2]).split('.')[0]+'.'+str(coordinates[2]).split('.')[1][:3]
    if 'e' in str(coordinates[2]):
        if int(str(coordinates[2]).split('e')[1]) < 0: z = '0.000'
    while len(x) <= 6: x = ' '+x
    while len(y) <= 7: y = ' '+y
    while len(z) <= 7: z = ' '+z
    while len(residue_number) <= 3: residue_number = ' '+residue_number
    while len(atom_number) <= 4: atom_number = ' '+atom_number
    while len(atom_name) <= 2: atom_name = atom_name+' '

    return 'ATOM  '+atom_number+'  '+atom_name+' '+residue_name+' '+chain+residue_number+'     '+x+y+z+'  0.00 00.00           '+atom_name[0]+'\n'

def range_find(coord, influence, boxside):

    if coord+influence > 0: top = int(coord+influence)+1
    else: top = int(coord+influence)
    
    if coord-influence > 0: bot = int(coord-influence)+1
    else: bot = int(coord-influence)

    if bot < -((boxside/2)-1) and top < -((boxside/2)-1): top = bot = -(boxside/2)
    elif bot < -((boxside/2)-1): bot = -((boxside/2)-1)
    if top > boxside/2 and bot > boxside/2: top = bot = (boxside/2)+1
    elif top > boxside/2: top = boxside/2

    indexbot = ((boxside/2)-1)+bot
    indextop = ((boxside/2)-1)+top
    
    return int(indextop), int(indexbot)


def sorted_insertion(value, code, sorted_list):

    for pos in range(len(sorted_list)):
        if value > sorted_list[pos][0]:
            sorted_list = sorted_list[:pos]+[[value, code]]+sorted_list[pos:]
            break
    if [value, code] not in sorted_list: sorted_list.append([value, code])

    return sorted_list

def h_clustering(points):
    hierarchy = []

    cluster = {}
    for n1 in range(len(points)):
        cluster[n1] = {'center':points[n1], 'components':[points[n1]]}

    distances = {}
    for n1 in range(len(points)):
        distances[n1] = []
        for n2 in range(len(points)):
            if n1 == n2: continue
            distance = distance_find(cluster[n1]['center'], cluster[n2]['center'])
            distances[n1] = sorted_insertion(distance, n2, distances[n1])

    min_distance = 10
    for key in distances:
        if min(min_distance, distances[key][-1][0]) == distances[key][-1][0]:
            min_distance = distances[key][-1][0]
            n2 = distances[key][-1][1]
            n1 = key

    while min_distance < 10:
 
        hierarchy.append([n1, n2])

        cluster[n1]['components'].extend(cluster[n2]['components'])
        accX, accY, accZ = [0,0,0]
        for coord in cluster[n1]['components']:
            accX += coord[0]
            accY += coord[1]
            accZ += coord[2]
        avgX = accX/len(cluster[n1]['components'])
        avgY = accY/len(cluster[n1]['components'])
        avzZ = accZ/len(cluster[n1]['components'])
        cluster[n1]['center'] = (avgX, avgY, avzZ)
        del distances[n1]
        del distances[n2]
        del cluster[n2]


        distances[n1] = []
        for key in distances:
            for distance in distances[key]: 
                if distance[1] == n2: distances[key].remove(distance)
            for distance in distances[key]: 
                if distance[1] == n1: distances[key].remove(distance)

        for key in distances:
            distance = distance_find(cluster[n1]['center'], cluster[key]['center'])
            if key == n1: continue
            distances[n1] = sorted_insertion(distance, key, distances[n1])
            distances[key] = sorted_insertion(distance, n1, distances[key])


        min_distance = 10
        for key in distances:
            if min(min_distance, distances[key][-1][0]) == distances[key][-1][0]:
                min_distance = distances[key][-1][0]
                n2 = distances[key][-1][1]
                n1 = key

    final_clusters = []
    for key in cluster: final_clusters.append(cluster[key]['components'])

    return final_clusters, hierarchy

def write_report(top_pred, pos_number, label):
    outfile = open('test_report', 'a')
    pos_top1x = 0
    pos_top5x = 0
    pos_top10x = 0
    for pos in top_pred[:pos_number]:
        if pos[1] == 1: pos_top1x += 1
    for pos in top_pred[:pos_number*5]:
        if pos[1] == 1: pos_top5x += 1
    for pos in top_pred[:pos_number*10]:
        if pos[1] == 1: pos_top10x += 1

    outfile.write(label+'('+str(pos_number)+' positives)\n')
    outfile.write('TPR top1x: '+str(pos_top1x/pos_number)+' ('+str(pos_top1x)+' in top'+str(pos_number)+' scores)\n')
    outfile.write('TPR top5x: '+str(pos_top5x/pos_number)+' ('+str(pos_top5x)+' in top'+str(pos_number*5)+' scores)\n')
    outfile.write('TPR top10x: '+str(pos_top10x/pos_number)+' ('+str(pos_top10x)+' in top'+str(pos_number*10)+' scores)\n\n')
    outfile.close()

###################################
##### STATISTICS CALCULATIONS #####
###################################

def metrics(cm):
##### given a confusion matrix in this format: {'PP':{'TP':0,'FP':0,},'PN':{'TN':0,'FN':0}}
##### returns a list with some common metrics

    POS=cm['PP']['TP']+cm['PN']['FN']
    NEG=cm['PN']['TN']+cm['PP']['FP']

    POSr=float(POS)/(POS+NEG)
    NEGr=float(NEG)/(POS+NEG)

    try: TPR=cm['PP']['TP']/POS
    except: TPR=0

    try: TNR=cm['PN']['TN']/NEG
    except: TNR=0

    try: PPV=cm['PP']['TP']/(cm['PP']['TP']+cm['PP']['FP'])
    except: PPV=0

    ACC=(cm['PP']['TP']+cm['PN']['TN'])/(cm['PP']['TP']+cm['PN']['TN']+cm['PP']['FP']+cm['PN']['FN'])

    MCCnum = (cm['PP']['TP']*cm['PN']['TN'])-(cm['PP']['FP']*cm['PN']['FN'])
    MCCden = (cm['PP']['TP']+cm['PP']['FP'])*POS*(cm['PN']['TN']+cm['PN']['FN'])*NEG
    try: MCC = MCCnum/math.sqrt(MCCden)
    except: MCC=0

    return [POS, NEG, TPR, 1-TNR, TNR, ACC, PPV, MCC]

def top_pred_ppv(label_list):
    test_cm = {'TP':0,'FP':0}

    for label in label_list:
        if label[1] == 1: test_cm['TP'] += 1
        else: test_cm['FP'] += 1
    
    return float(test_cm['TP'])/(float(test_cm['TP'])+float(test_cm['FP']))

def fisher_test(tab, top):

    ab = math.factorial(tab['t']['p']+tab['t']['n'])
    cd = math.factorial(tab['nt']['p']+tab['nt']['n'])
    ac = math.factorial(tab['t']['p']+tab['nt']['p'])
    bd = math.factorial(tab['t']['n']+tab['nt']['n'])
    fnum = ab*cd*ac*bd

    abcd =  tab['t']['p']+tab['t']['n']+tab['nt']['p']+tab['nt']['n']
    fden = math.factorial(tab['t']['p'])*math.factorial(tab['t']['n'])*math.factorial(tab['nt']['p'])*math.factorial(tab['nt']['n'])*math.factorial(abcd)

    f = fnum/fden

    if tab['nt']['p'] > 0 and tab['t']['p'] < top: 
        more_extreme_table = copy.deepcopy(tab)
        more_extreme_table['t']['p'] += 1
        more_extreme_table['t']['n'] -= 1
        more_extreme_table['nt']['p'] -= 1
        more_extreme_table['nt']['n'] += 1
        return f + fisher_test(more_extreme_table, top)

    else: return f 

#################################
##### KERAS LAYER FUNCTIONS #####
#################################

def dense (out, act, reg, drp, bn, inl):

    a = Dense(out, 
              activation=act, 
              use_bias='True', 
              kernel_regularizer=l2(reg))(inl)

    if bn: a = BatchNormalization(axis=-1)(a)

    a = Dropout(drp)(a)

    return a

def conv1D(out, ks, act, reg, drp, bn, inl):

    a = Conv1D(out, 
               kernel_size=(ks), 
               activation=act, 
               strides=1, 
               padding='valid', 
               use_bias='True', 
               kernel_regularizer=l2(reg))(inl)

    if bn: a = BatchNormalization(axis=-1)(a)

    a = Dropout(drp)(a)

    return a

def conv2D(out, ks, act, dr, reg, drp, bn, inl):

    a = Conv2D(out, 
               kernel_size=ks, 
               dilation_rate=dr,
               activation=act, 
               strides=1, 
               padding='valid', 
               use_bias='True',
               kernel_regularizer=l2(reg), 
               data_format='channels_last')(inl)

    if bn: a = BatchNormalization(axis=-1)(a)

    a = Dropout(drp)(a)

    return a

def conv3D(out, ks, act, dr, reg, drp, bn, inl):

    a = Conv3D(out, 
               kernel_size=ks,
               dilation_rate=dr,
               activation=act, 
               strides=1, 
               padding='valid', 
               use_bias='True',
               kernel_regularizer=l2(reg),  
               data_format='channels_last')(inl)

    if bn: a = BatchNormalization(axis=-1)(a)

    a = Dropout(drp)(a)

    return a

def bi_LSTM(out, reg, drp, seq, bn, inl):

    a = Bidirectional(CuDNNLSTM(out, 
                                return_sequences=seq, 
                                kernel_regularizer=l2(reg),
                                recurrent_regularizer=l2(reg), 
                                activity_regularizer=l2(reg)),
                      merge_mode='concat')(inl)

    if bn: a = BatchNormalization(axis=-1)(a)

    a = Dropout(drp)(a)

    return a

def bi_GRU(out, reg, drp, seq, bn, inl):

    a = Bidirectional(CuDNNGRU(out, 
                               return_sequences=seq, 
                               kernel_regularizer=l2(reg),
                               recurrent_regularizer=l2(reg), 
                               activity_regularizer=l2(reg)),
                      merge_mode='concat')(inl)

    if bn: a = BatchNormalization(axis=-1)(a)

    a = Dropout(drp)(a)

    return a

README - FILE 
Naccess - accessibility calculations
------------------------------------
Simon Hubbard, 	Biomolecular Structure and Modelling Unit, 
		University College, 
		Gower Street, 
		London WC1E 6BT
                UK. 
(present address):
                Department of Biomolecular Sciences,
                UMIST, PO Box 88
                Manchester M60 1QD
                UK.

******************************************************************************
* PLEASE NOTE - THERE IS A CONFIDENTIALITY AGREEMENT ATTACHED TO THE BOTTOM  *
* OF THIS README FILE. IF YOU HAVEN'T DONE SO ALREADY, PLEASE PRINT IT, SIGN *
* IT, AND SEND IT TO SIMON HUBBARD, AT THE PRESENT ADDRESS GIVEN ABOVE.      *
* THANKS VERY MUCH FOR YOUR COOPERATION                                      *
*****************************************************************************

VERSION 2.1.1 S.J.Hubbard, June 1996.

Briefly, the naccess program calculates the atomic accessible surface defined
by rolling a probe of given size around a van der Waals surface. This program
is an implimentation of the method of Lee and Richards (1971) J.Mol.Biol.55,
379-400. which does just that. The program is dimensioned for up to 20000
atoms, and allows the variation of the probe size and atomic radii by the user.
The program is written in (fairly) standard FORTRAN 77 and should compile
on most UNIX platforms. It outputs 3 files: 

1) An atomic accessibility file (.asa file) containing the calculated
accessible surface for each atom in a PDB file, as well as the assigned van der
Waal radii.

2) A residue accessibility (.rsa) file containing summed atomic accessible
surface areas over each protein or nucleic acid residue, as well as the
relative accessibility of each residue calculated as the %accessiblity compared
to the accessibility of that residue type in an extended ALA-x-ALA
tripeptide (for amino acids). See Hubbard, Campbell & Thornton (1991) 
J.Mol.Biol.220,507-530. You can prevent this file from being calculated 
if you don't need it.

3) A log file (.log) containing information concerning the calculation.


*********************************UPDATE 13/6/96:********************************
**NEW** Now handles nucleic acids. 5 standard (A,C,G,T,U) included in vdw.radii
**NEW** New .rsa output format, 7 columns instead of 5, with 2 additional fields
        for All-polar and all-non-polar atoms
**NEW** CHAIN records containing summations over each individual chain
**NEW** -q option prints the options list
********************************************************************************

installation
------------
If you are reading this file you may have already done it ! You should
have the following files in your installation directory:

README          accall.pars     naccess.scr*    vdw.radii
accall.f        install.scr     standard.data

next just type "csh install.scr", and the program should be compiled
for you, and then final script that runs everything should be "naccess"

usage
-----
The program is run by a script called "naccess". At installation, the script
is aliased to the command "naccess", but you should create your own alias, and
include it in your .cshrc or some other script if you want to use the program
regularly. Something of the form:

alias naccess '/home/username/naccess_directory/naccess'

To use the program, simply type "naccess", and the name of a valid PDB file,
including the full path where appropriate. There are a number of other
parameters that can be supplied, but in its simplest form, it might be
something like:

naccess /data/pdb/1crn.brk

and would produce the files:

1crn.asa 1crn.log 1crn.rsa

The full usage is:
naccess pdb_file -p probe_size -r vdw_radii_file -s std_data_file -z zslice -[hwyfaclqb]"
NOTE: use multiple options separately, ie: naccess test.pdb -p 1.20 -h -f -a

The default probe size is 1.40 Angstroms. To use a probe of size 1.20 A type: 

naccess /data/pdb/1crn.brk -p 1.2

If you don't want to use the default van der Waal radii, taken from Chothia
(1976) J.Mol.Biol.105,1-14, then you copy the file "vdw.radii" to another file
and edit the values appropraitely. Then you can type:

naccess /data/pdb/1crn.brk -r my.radii

If you supply a file which doesn't exist, the program will default to the
vdw.radii file. It looks for it first in the current directory, then the
naccess executable directory, and then gives up !!

The Lee and Richards method works by taking thin Z-slices through the molecule
and calculating the exposed arc lengths for each atom in each slice, and then
summing the arc lengths to the final area over all z-values. Hence, the zslice
parameter controls accuracy and also speed of calculation. The default value is
0.05 A, but can be changed. For rough, quick calculations, a value of 0.1 might
be better. You can type:

naccess /data/pdb/1crn.brk -z 0.1

As a rough estimate, the program takes about 8 secs (real time) for crambin
(1crn) using the default parameters on an IRIS R4000 Indigo workstation. Thats
less than 2 cpu seconds.

By default, the program ignores HETATM records, hydrogens, and waters. If you
want these to be considered in the calculation supply a parameter of the form 
-h, -w and/or -y respectively. You can type:

naccess /data/pdb/5pti.brk -y 

  or

naccess /data/pdb/4pti.brk -w 

  or

naccess /data/pdb/4hhb.brk -h

to see the effect of these parameters on the output files.

Naccess now handles nucleic acids explicitly and the standard 5 (A,C,G,T,U) are
included in the vdw.radii file supplied with the distribution. Add your own by
editting this file and maintaining the format.

If you require atomic *contact* areas rather than accessible areas - this
is not the path traced by the probe centroid but the parts of the vdw surface
that the probe actually touches - then use the -c option:

naccess /data/pdb/4pti.brk -c 

NOTE: The relative accessibilities (see .rsa files below) will be incorrect if
you use this option. You can relcalulate them as contact areas if you wish 
with standard ala-x-ala tripeptides, and create a new "standard.data" type 
file to rectify this. I have the data for this if required - just send email.

Apart from HEME groups, the van der Waals radii are not explicitly defined in
the default vdw.radii file. You should add your own in if there are other
HETATM groups you are explicitly interested in. Otherwise the program makes
crude guesses at the respective radii !!

The .rsa residue accessibility file is created using the file "standard.data"
to calculate the percentage accessibilities. If these are unsatisfactory, edit
them. Again, the program looks in the current directory for this file, and then
the naccess executable directory. You can supply an alternative file using the
-s option, which works like the -r option for van der Waals radii.

NOTE: The default .rsa output format has changed. The old format (plus 2 new
fields) can be obtained using the -l (l for long) format. See below.

By default, alpha carbon atoms are considered to be part of the amino acid side 
chain so that glycines possess a relative side chain accessibility. To switch 
this off, use the -b option. WARNING: Please note that when using this option, 
the sidechain and mainchain %accessibilities for all residues will be wrong.

multiple chains
---------------

Naccess calcualtes the accessibilities of the whole molecule submitted in a 
PDB file. Hence, it is down to the user to control the input if they want to
calculate the accessibility of a single chain from a multi-PDB file. This is a 
conscious decision that allows total flexibility, and prevents the program
from interfering with the co-ordinates, perhaps unexpectedly. Hence, at present,
if users wish to calculate the accessible area buried between 2 chains, they
will have to run 3 separate calculations. Once for chain A, once for chain B
and a third time for the AB complex. Then they will have to subtract the
AB surface from the A + B surfaces to obtain the buried surface in the complex.
Hopefully in the next version of Naccess this will be automated.

example output files
--------------------

The first 2 residues from an example .asa file are shown below. The output
format is PDB, with B-factors and occupancies removed, then atomic accessiblity
in square Angstroms, followed by the assigned van der Waal radius. If you
want to keep the occupancies and B-factors in, then use the -f option (full)
which gives a slightly different (and larger) output format, but is compatible
with the hydrogen bond calculation program HBPLUS from Ian McDonald.

ATOM      1  N   THR     1      17.047  14.099   3.625  22.279  1.65
ATOM      2  CA  THR     1      16.967  12.784   4.338  13.902  1.87
ATOM      3  C   THR     1      15.685  12.755   5.133   0.000  1.76
ATOM      4  O   THR     1      15.268  13.825   5.594   0.000  1.40
ATOM      5  CB  THR     1      18.170  12.703   5.337   0.098  1.87
ATOM      6  OG1 THR     1      19.334  12.829   4.463  17.632  1.40
ATOM      7  CG2 THR     1      18.150  11.546   6.304  20.662  1.87
ATOM      8  N   THR     2      15.115  11.555   5.265   3.568  1.65
ATOM      9  CA  THR     2      13.856  11.469   6.066   0.000  1.87
ATOM     10  C   THR     2      14.164  10.785   7.379   0.000  1.76
ATOM     11  O   THR     2      14.993   9.862   7.443   5.072  1.40
ATOM     12  CB  THR     2      12.732  10.711   5.261   0.010  1.87
ATOM     13  OG1 THR     2      13.308   9.439   4.926  10.623  1.40
ATOM     14  CG2 THR     2      12.484  11.442   3.895   2.739  1.87

The beginning and end of an example .rsa file are shown below. The data is
summed over residues, and split into 5 classes. Total (all atoms), Non Polar
Sidechain (all non-oxygens and non-nitrogens in the sidechain), Polar Sidechain
(all oxygens and nitrogens in the sidechain), total sidechain, and
mainchain. For our purposes, alpha carbons are classed as sidechain atoms, so
that glycine can have a sidechain accessibility. They are therefore not
included in the mainchain. For each class, two values are given, an absolute
(ABS) and relative (REL) accessibility. The absolute value is the simple sum,
whilst the REL value is the % relative accessiblity. Absolute sums over the
whole chain are also given.

To avoid calculating and producing such a summary file, supply the program
with the -a option, which speeds things up very slightly.

REM  Relative accessibilites read from external file "standard.data"
REM  File of summed (Sum) and % (per.) accessibilities for 
REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
RES ARG     1   134.65  56.4  78.55  39.0  56.10 149.6  13.44  17.3 121.22  75.3
RES PRO     2    50.10  36.8  49.91  41.6   0.19   1.2  49.91  41.3   0.19   1.2
RES ASP     3   124.60  88.8 112.69 109.7  11.91  31.6  43.08  87.5  81.52  89.4
RES PHE     4    37.07  18.6  35.12  21.4   1.95   5.5  35.12  21.3   1.95   5.7
RES CYS     5     0.16   0.1   0.00   0.0   0.16   0.4   0.00   0.0   0.16   0.5
RES LEU     6   108.00  60.5  83.45  59.1  24.55  65.5  85.45  60.0  22.55  62.1
RES GLU     7    52.21  30.3  49.52  36.8   2.69   7.2  27.45  45.5  24.75  22.1
RES PRO     8   102.59  75.4 100.85  84.1   1.74  10.7 100.85  83.4   1.74  11.5
RES PRO     9    53.39  39.2  34.67  28.9  18.72 115.3  34.87  28.8  18.52 121.9
RES TYR    10    79.16  37.2  76.48  43.1   2.68   7.6  48.91  35.8  30.25  39.7
 .   .      .      .      .     .      .     .      .     .      .     .      .
 .   .      .      .      .     .      .     .      .     .      .     .      .
 .   .      .      .      .     .      .     .      .     .      .     .      .
RES MET    52    72.80  37.5  65.18  41.6   7.61  20.3  65.92  41.8   6.87  18.9
RES ARG    53   178.12  74.6 154.31  76.7  23.80  63.5  80.72 103.8  97.40  60.5
RES THR    54    49.34  35.4  42.01  41.3   7.33  19.5  23.21  30.7  26.13  41.1
RES CYS    55     0.00   0.0   0.00   0.0   0.00   0.0   0.00   0.0   0.00   0.0
RES GLY    56    27.18  33.9   3.02   9.4  24.15  50.6   6.18  16.5  21.00  49.3
RES GLY    57    34.25  42.8   5.23  16.2  29.02  60.7   9.85  26.2  24.40  57.3
RES ALA    58   167.77 155.4  82.34 118.6  85.44 221.7  89.05 124.8  78.72 215.2
END  Absolute sums over single chains surface 
CHAIN  1 _     3971.5       3280.0        691.4       2209.7       1761.8
END  Absolute sums over all chains 
TOTAL          3971.5       3280.0        691.4       2209.7       1761.8

An example .log file follows

 ACCALL - Accessibility calculations
 PDB FILE INPUT /data/pdb/1crn.brk
 PROBE SIZE       1.40
 Z-SLICE WIDTH   0.050
 VDW RADII FILE vdw.radii
 EXCL HETATOMS
 EXCL HYDROGENS
 EXCL WATERS
 READVDW  25 residues input
 NON-STANDARD atom. OXT in residue> ASN    46 
 GUESSED vdw of  OXT in ASN    46  =  1.40
 ADDED VDW RADII
 CHAINS       1
 RESIDUES    46
 ATOMS      327
 SOLVA: PROGRAM ENDS CORRECTLY
 CALCULATED ATOMIC ACCESSIBILITES
 RELATIVE (STANDARD) ACCESSIBILITIES READFOR  23 AMINO ACIDS
 SUMMED ACCESSIBILITIES OVER RESIDUES

support
-------
All queries and questions to "sjh@sjh.bi.umist.ac.uk". I hope it all works.

-----------------------------------------------------------------------
-----------------------------------------------------------------------


CONFIDENTIALITY AGREEMENT
=========================

Correspondence to:
Dr. Simon Hubbard
Department of Biomolecular Sciences,
UMIST, PO Box 88,
Manchester M60 1QD, UK

Email (Internet): sjh@sjh.bi.umist.ac.uk


	           Naccess - Accessibility calculations
                   ------------------------------------          

			CONFIDENTIALITY AGREEMENT
			-------------------------


In regard to the NACCESS programs, specified in Appendix 1 herewith (the
Software) supplied to us, the copyright and other intellectual property
rights to which belong to the authors, we

    __________________________________________________________________

undertake to the authors that we shall be bound by the following terms and
conditions:-

1. We will receive the Software and any related documentation in confidence
and will not use the same except for the purpose of the department's own 
research. The Software will be used only by such of our officers or
employees to whom it must reasonably be communicated to enable us to
undertake our research and who agree to be bound by the same confidence.
The department shall procure and enforce such agreement from its staff for
the benefit of the authors.

2. The publication of research using the Software should reference Hubbard,S.J.
& Thornton, J.M. (1993), 'NACCESS', Computer Program, Department of 
Biochemistry and Molecular Biology, University College London." or successor 
references as defined by the authors.

3. Research shall take place solely at the department's premises at

    __________________________________________________________________

4. All forms of the Software will be kept in a reasonably secure place to
prevent unauthorised access.

5. Each copy of the Software or, if not practicable then, any package
associated therewith shall be suitably marked (and such marking maintained)
with the following copyright notice: "Copyright 1992-3 Simon Hubbard and 
Janet M Thornton All Rights Reserved".

6. The Software may be modified but any changes made shall be made
available to the authors.

7. The Software shall be used exclusively for academic teaching and
research. The Software will not be used for any commercial research or
research associated with an industrial company.

8. The confidentiality obligation in paragraph one shall not apply:

   (i)  to information and data known to the department at the time of
	receipt hereunder (as evidenced by its written records);

  (ii)	to information and data which was at the time of receipt in the 
	public domain or thereafter becomes so through no wrongful act of
	the department;

 (iii)	to information and data which the department receives from a third
	party not in breach of any obligation of confidentiality owed to
	the authors.



Please sign this Undertaking and return a copy of it to indicate that you 
have read, understood and accepted the above terms.



		      For and on behalf of _____________________________

		      _________________________________________________
		     
		      ..................................................

		      Dated ............................................

-------------------------------------------------------------
Appendix 1 -  Files supplied as part of NACCESS

Source code:
accall.f

Script files:
install.scr    
naccess.scr

Data files:
standard.data
accall.pars         
vdw.radii

Documentation:
README        

README - END


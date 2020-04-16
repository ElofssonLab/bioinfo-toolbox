###################### REQUIREMENTS ##################################
### To do before the usage:
### - download tr_Rosetta (git clone https://github.com/gjoni/trRosetta) [for these instructions is called trRosetta_1]
### - download pre-trained network and untar (wget https://files.ipd.uw.edu/pub/trRosetta/model2019_07.tar.bz2 ; tar xf model2019_07.tar.bz2)
### - download second part of trRosetta (server), via registration (http://yanglab.nankai.edu.cn/trRosetta/) [for these instructions is called trRosetta_2]
###
### - build the singularity container:
###       download Pyrosetta zip Ubuntu 18.04 LTS (64-bit) and Python-3.6.Release (http://www.pyrosetta.org/dow)
###       singularity build --sandbox --fakeroot directoryname tf1.def
###       singularity shell --writable directoryname
###################### USAGE ##################################
### generate an MSA (in a3m format) for your protein sequence
### be sure that the singularity container contains both the tr_Rosetta packages and fasta + MSA of your protein sequence
### 1) python trRosetta_1/network/predict.py -m ./model2019_07 fastaseq.a3m fastaseq.npz
### 2) python trRosetta_2/trRosetta.py fastaseq.npz fastaseq.fasta fastaseq.pdb
### it is suggested running step two for multiple times to generate multiple models and select the top models based on the energy scores, which are available at the end of the model's pdb file.
### More details about trRosetta can be found from the following paper:
###J Yang et al, Improved protein structure prediction using predicted inter-residue orientations, PNAS (2020).

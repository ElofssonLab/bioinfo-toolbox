# Short instructions on how to use trRosetta for modelling

## Installation of trRosetta, PyRosetta and Rosetta 


*  download tr_Rosetta
  * git clone https://github.com/gjoni/trRosetta


* 
### - download pre-trained network and untar (wget
###   https://files.ipd.uw.edu/pub/trRosetta/model2019_07.tar.bz2 ;
###   tar xf model2019_07.tar.bz2)


### - download second part of trRosetta (server), via registration
###   (http://yanglab.nankai.edu.cn/trRosetta/) [for these
###   instructions is called trRosetta_2]


###
### - build the singularity container:
###       download Pyrosetta zip Ubuntu 18.04 LTS (64-bit) and
###       Python-3.6.Release (http://www.pyrosetta.org/dow)
###       Note - make sure you get the right Python version
###       singularity build --sandbox --fakeroot directoryname sing1.def
###       singularity shell --writable directoryname
###       or 
###       sudo singularity build  sing1.def

### Command used by Arne (as he was runnin gout of memory). To build GPU enabled image
### sudo su -c "SINGULARITY_TMPDIR=/home/arnee/Downloads/tmp/  SINGULARITY_DISABLE_CACHE=true SINGULARITY_CACHE=/home/arnee/Downloads/tmp/ singularity build /home/arnee/singularity-images/trRosetta-gpu.simg sing1-gpu.def"


###################### USAGE ##################################
### generate an MSA (in a3m format) for your protein sequence
### be sure that the singularity container contains and the tr_Rosetta
### packages and fasta + MSA of your protein sequence are accesible from the
### container

### 1) Generate the distance map
### > python trRosetta_1/network/predict.py -m ./model2019_07 fastaseq.a3m fastaseq.npz

### 2) Generare one model
### > python trRosetta_2/trRosetta.py fastaseq.npz fastaseq.fasta fastaseq.pdb

###    it is suggested running step two for multiple times to generate
###    multiple models and select the top models based on the energy
###    scores, which are available at the end of the model's pdb file.

### More details about trRosetta can be found from the following
### paper:J Yang et al, Improved protein structure prediction using
### predicted inter-residue orientations, PNAS (2020).


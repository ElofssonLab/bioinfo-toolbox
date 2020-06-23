# Short instructions on how to use trRosetta for modelling

## Installation of trRosetta, PyRosetta and Rosetta


*  download tr_Rosetta
```bash
git clone https://github.com/gjoni/trRosetta
```
*  download pre-trained network and untar
```bash
wget https://files.ipd.uw.edu/pub/trRosetta/model2019_07.tar.bz2
tar xf model2019_07.tar.bz2)
```
* download second part of trRosetta (server)
  * registration necessary at http://yanglab.nankai.edu.cn/trRosetta/)

* build the singularity container (Also available on kebnekaise and on Arnes machine)
  * Two definitions available here
    * Default sing1.def (with ipython) and sing1-gpu.def (with Cuda)
    * singularity build directoryname sing1.def  (or use fakeroot/sannbox if you prefer)
      * Command used by Arne (as he was running out of memory). To build GPU enabled image
       ```bash
       sudo su -c "SINGULARITY_TMPDIR=/home/arnee/Downloads/tmp/  SINGULARITY_DISABLE_CACHE=true SINGULARITY_CACHE=/home/arnee/Downloads/tmp/ singularity build /home/arnee/singularity-images/trRosetta-gpu.simg sing1-gpu.def"
       ```

* download Pyrosetta zip file for Ubuntu 18.04 LTS (64-bit) and Python-3.6  from
  http://www.pyrosetta.org/dow
  * Note:  make sure you get the right Python version


## List of files with short explanations (note things are not ready for publication)
* Lets assume we whan to "dock X and Y" (fasta files X.fa, Y.fa)
* To prepare the MSA*
```bash
  jackhmmer -N 3 -o X.jackhmmer -A X.sto X.fa database.fasta
  jackhmmer -N 3 -o Y.jackhmmer -A Y.sto Y.fa database.fasta
  python3 mergeSTObyGenome.py -i X.sto -j Y.sto > X-Y.a3m
  python3 a3mToTrimmed.py  X-Y.a3m > X-Y.trimmed
```
\* "Trim an a3m file - probably not needed

* Run trRosetta distance predictionsx*
```bash
  python3 network/predict.py  -m ~/model2019_07 X-Y.trimmed X-Y.npz
```
\* it might also be worth to try to run Y-X.a3m (reversed order)

* Modify the NPX file (to delete the insertion
 ```bash 
  python3 del_sepseq_npz.py -i X-Y.npz -s X-Y.fasta -o X-Y-nosep.npz
 ```
   \- Modify the NPZ file by deleting the npz file.

* To run pyRosetta
```bash
  python3 trRosetta-2chain.py X-Y-nodel.npx X.fa Y.fa X-Y.pdb
```
   \- The script used to fold two chains.
   
  * We need (perhaps) to limit the constraints that attracts the two sequences
    to only the ones that have prob>0.5
    - not done yet.

* Utils needed by pyrosetta
  * arguments.py
    - one extra argument to handle two fasta input files)
  * data/
    - Files with potentials
  * utils_ros.py
    - New subrouting necessary

* Some useful programs
  * npz_to_rr.py
    - Convert and trRosettta predicton file to a CASP contact file
  * sto_to_msa.py
    - Convert stockholm to a3m (fasta format
  * tr_calculate_contact_dist.py
    - makes a plot of the distance map (with poly G shaded)
  * start.py
    - The first part of trRosetta-2chain.py
    - useful for ipython



# Some extra info
*   it is suggested running step two multiple times to generate
   several models and select the top models based on the energy
    scores, which are available at the end of the model's pdb file.

* More details about trRosetta can be found from the following
 paper:J Yang et al, [Improved protein structure prediction using
 predicted inter-residue orientations](https://www.pnas.org/content/117/3/1496), PNAS (2020).



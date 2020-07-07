# Short instructions on how to use trRosetta for modeling

## Installation of PyRosetta and trRosetta

*  download tr_Rosetta
```bash
git clone https://github.com/gjoni/trRosetta
```
*  download pre-trained network and untar
```bash
wget https://files.ipd.uw.edu/pub/trRosetta/model2019_07.tar.bz2
tar xf model2019_07.tar.bz2
```
* download second part of trRosetta (server)
  * registration necessary at http://yanglab.nankai.edu.cn/trRosetta/)

* download Pyrosetta zip file for Ubuntu 18.04 LTS (64-bit) and Python-3.6 from
  http://www.pyrosetta.org/download
  Note:  make sure you get the right Python version

* --- for Snakemake pipeline --- 
  * download CNS: https://www.dropbox.com/s/748qg15ywq2vfv8/cns_solve_1.3_all_intel-mac_linux.tar.gz?dl=0
  * make sure to copy cns_solve_env.sh from this folder to the folder where you launch the singularity build command

* build the singularity container (Also available on kebnekaise and on Arnes machine)
  * Three definition files available here
    * Default sing1.def (with ipython), sing1-gpu.def (with Cuda) and sing1-snake.def (with all the Snakemake executables)
  * To build the singularity image:
  ```bash
  singularity build --fakeroot image_name sing1.def  (also use --sandbox if you want to be able to modify the container later)
  ``` 
  * To build GPU enabled image in case of out of memory error:
  ```bash 
  sudo su -c "SINGULARITY_TMPDIR=/home/arnee/Downloads/tmp/  SINGULARITY_DISABLE_CACHE=true SINGULARITY_CACHE=/home/arnee/Downloads/tmp/ singularity build /home/arnee/singularity-images/trRosetta-gpu.simg sing1-gpu.def"
  ```

## List of files with short explanations (note things are not ready for publication)
* Lets assume we want to "dock X and Y" (fasta files X.fa, Y.fa)
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

* Generate a "fasta" file that contains the sequence used predicting the distances. (make sure you get the whole sequence)
 ```bash 
  head -2 X-Y.trimmed > X-Y.fasta
 ```
   \- Modify the NPZ file by deleting the npz file.
* Modify the NPZ file (to delete the insertion). 
 ```bash 
  python3 del_sepseq_npz.py -i X-Y.npz -s X-Y.fasta -o X-Y-nosep.npz
  ```
   \- In addition one can merge the files from two individual predictions and
 the merged one.
* In addition merge the files from two individual predictions and
 the merged one. 
 ```bash 
  python3 merge_npz.py -i X.npz -j Y.npz  -k X-Y-nosep.npz -o X-Y-merged.npz
  ```
   \- Modify the NPZ file by deleting the npz file.

* To run pyRosetta
```bash
  python3 trRosetta-2chain.py X-Y-model.npz X.fa Y.fa X-Y.pdb
```
	* Some parameters 
	  * -minprob the minimum probability of interactions to be included.
	  * -intradist  (distance for the "flat" potential, i.e. where attraction should start. Default 15 (as when it is closer the real potential can start to act)
	  * -intrasd  Standard deviation of the intra-chain potential (default 20)
	  * -allintra  (adds an even weaker constraint on all intra pairs)

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
    - New subroutine necessary

* Some useful programs
  * npz_to_rr.py
    - Convert and trRosetta prediction file to a CASP contact file
  * sto_to_msa.py
    - Convert stockholm to a3m (fasta format
  * analyze_npz.py
    - makes a plot of the distance map (with poly G shaded)
    - Calcualates statistics of the contact map -useful for scoring.
  * start.py
    - The first part of trRosetta-2chain.py
    - useful for ipython


# Snakemake  (Copy to your working directory)
* There is a snakemake file that makes it possible to run the while pipeline in one step.
  * Direcotory structure:
    - seq/X.fasta   - single sequence fasta files Necessary input
    - msa  - Jackhmmer output
    - dimer/ - merged MSA
    - distpred/ - NPZ file from trRosetta
    - results/ - some analyzis and plots of distance map
    - model/ - output PDB file (currently only one)
  * To run
   ```bash 
   	  snakemake model/X-Y.pdb   (model/Y-X.pdb will use the other sequential order)
   ```
  * If you want to have the analysis also run
   ```bash 
   	  snakemake results/X-Y.png 
   ```
  
  

# Some extra info
*   it is suggested running step two multiple times to generate
   several models and select the top models based on the energy
    scores, which are available at the end of the model's pdb file.

* More details about trRosetta can be found from the following
 paper:J Yang et al, [Improved protein structure prediction using
 predicted inter-residue orientations](https://www.pnas.org/content/117/3/1496), PNAS (2020).



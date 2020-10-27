# there are two main scripts to run: blast_rbh.sh and jackhmmer_rbh.sh. 

**blast_rbh.py**:

  - modify the parameters in *param setting* section. This section includes the setting of: project folder, bioinfo-toolbox, ecoli fasta file dir(ecoli fasta blast database dir if using PSIBLAST), ecoli ID list txt file, and the paths of saving different outputs.
  - PSIBLAST method can only search homologs against reference proteomes(8841 proteomes). 
  - PSIBLAST command is : 3 iterations, evalue 0.01. 


**jackhmmer_rbh.py**:
  - modify the parameters in *param setting* section. It's similar to *blast_rbh*.
  - Jackhmmer method can search homologs against reference and nonreference proteomes(32979 proteomes).
  - Jackhmmer command is : 3 iterations, evalue 0.01.





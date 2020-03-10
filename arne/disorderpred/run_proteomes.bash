#!/bin/bash -x

for i in proteomes/*_DNA.fasta
do	 

	 /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -rnafile -pro -f proteomes/UP000001037_694429_DNA.fasta -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_54.ann

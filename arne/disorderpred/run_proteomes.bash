#!/bin/bash -x

for i in proteomes/*_DNA.fasta
do	 
    j=`basename $i _DNA.fasta`
    # Protein
#    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py  -pro -f proteomes/$j.fasta -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_54.ann > proteomes/$j.fasta.model_pro
#    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py  -pro -f proteomes/$j.fasta -m models/model_conv1d_100-test1234.sorted.GCgenomic-10-0.001_1_pro_70.ann > proteomes/$j.fasta.model_pro_conv1d
#    
#    # Protein + GC
#    
#    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -rnafile -pro -gc -f $i -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_16.ann >$i.model_pro_GC
#    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -rnafile -pro -gc -f $i -m models/model_conv1d_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_78.ann >$i.model_conv1d_pro_GC

    # Protein + GCgenomiv

    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -pro -gcgenomic -f proteomes/$j.fasta -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_95.ann >  proteomes/$j.fasta.model_pro_GCgenomic    
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -pro -gcgenomic -f proteomes/$j.fasta -m models/model_conv1d_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_9.ann >  proteomes/$j.fasta.model_conv1d_pro_GCgenomic    

done

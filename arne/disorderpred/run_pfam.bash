#!/bin/bash -x

for dir in pfam/*/.
 
do
    for i in $dir/*mingap*.fa
    do
	# Protein
	/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py  -pro -f $i -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_54.ann > $i.fasta.model_pro
	/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py  -pro -f $i -m models/model_conv1d_100-test1234.sorted.GCgenomic-10-0.001_1_pro_70.ann > $i.model_pro_conv1d
	
	# Protein + GC
	#    
	#    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -rnafile -pro -gc -f $i -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_16.ann >$i.model_pro_GC
	#    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -rnafile -pro -gc -f $i -m models/model_conv1d_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_78.ann >$i.model_conv1d_pro_GC
	
	# Protein + GCgenomiv
	
	/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -pro -gcgenomic -f $i -m models/model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_95.ann >$i.model_pro_GCgenomic    
	/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_pred.py -pro -gcgenomic -f $i -m models/model_conv1d_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_9.ann >$i.model_conv1d_pro_GCgenomic    
    done
done

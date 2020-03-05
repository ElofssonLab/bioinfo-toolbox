#!/bin/bash -x


#python3 bin/mobiformatter.py -pro -f mobidata_K.pickly
#python3 bin/mobiformatter.py -pro -gc -f mobidata_K.pickly
#python3 bin/mobiformatter.py -pro -gcgenomic -f mobidata_K.pickly
#python3 bin/mobiformatter.py -pro -gc -gcgenomic -f mobidata_K.pickly

#python3 bin/mobiformatter.py -rna -f mobidata_K.pickly
#python3 bin/mobiformatter.py -rna -gc -f mobidata_K.pickly
#python3 bin/mobiformatter.py -rna -gcgenomic -f mobidata_K.pickly
#python3 bin/mobiformatter.py -rna -gc -gcgenomic -f mobidata_K.pickly

#For simplicity use test1-3 for training and test4-5 for testing  (Perhaps a too small set ?)



/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1-3 -v DIScv/test4-5 -d formatted_data_GC.h5py -f pro -gc
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1-3 -v DIScv/test4-5 -d formatted_data_GC.h5py -f rna -gc
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1-3 -v DIScv/test4-5 -d formatted_data_GC.h5py -f pro 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1-3 -v DIScv/test4-5 -d formatted_data_GC.h5py -f rna 

#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f pro -m logs/model_200-1-0.001_1_pro
#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f pro -m logs/model_200-1-0.001_1_rna


# Testing.

for i in models/model*rna*noGC*
do
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f rna -m $i
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/formatted_list -d ./ -f rna -m $i
done

for i in models/model*pro*noGC*
do
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f pro -m $i
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/formatted_list -d ./ -f pro -m $i
done


for i in models/model*rna*_GC*
do
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f rna -gc -m $i
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/formatted_list -d ./ -f rna -gc -m $i
done

for i in models/model*pro*_GC*
do
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f pro -gc -m $i
    /opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/formatted_list -d ./ -f pro -gc -m $i
done



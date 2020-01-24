#!/bin/bash -x


#python3 bin/mobiformatter.py

#For simplicity use test1-3 for training and test4-5 for testing  (Perhaps a too small set ?)



/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1-3 -v DIScv/test4-5 -d formatted_data_GC.h5py -f pro -gc
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1-3 -v DIScv/test4-5 -d formatted_data_GC.h5py -f rna -gc

#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f pro -m logs/model_200-1-0.001_1_pro
#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test.py  -t DIScv/test4-5 -d ./ -f pro -m logs/model_200-1-0.001_1_rna



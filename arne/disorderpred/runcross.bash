#!/bin/bash -x


/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1234 -v DIScv/test5 -d formatted_data_GC.h5py -f pro -gc -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1235 -v DIScv/test4 -d formatted_data_GC.h5py -f pro -gc -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1245 -v DIScv/test3 -d formatted_data_GC.h5py -f pro -gc  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1345 -v DIScv/test2 -d formatted_data_GC.h5py -f pro -gc -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test2345 -v DIScv/test1 -d formatted_data_GC.h5py -f pro -gc  -bs 10


/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1234 -v DIScv/test5 -d formatted_data_GC.h5py -f pro  -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1235 -v DIScv/test4 -d formatted_data_GC.h5py -f pro  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1245 -v DIScv/test3 -d formatted_data_GC.h5py -f pro   -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1345 -v DIScv/test2 -d formatted_data_GC.h5py -f pro  -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test2345 -v DIScv/test1 -d formatted_data_GC.h5py -f pro   -bs 10


/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1234 -v DIScv/test5 -d formatted_data_GC.h5py -f rna -gc -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1235 -v DIScv/test4 -d formatted_data_GC.h5py -f rna -gc -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1245 -v DIScv/test3 -d formatted_data_GC.h5py -f rna -gc  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1345 -v DIScv/test2 -d formatted_data_GC.h5py -f rna -gc -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test2345 -v DIScv/test1 -d formatted_data_GC.h5py -f rna -gc  -bs 10


/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1234 -v DIScv/test5 -d formatted_data_GC.h5py -f rna  -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1235 -v DIScv/test4 -d formatted_data_GC.h5py -f rna  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1245 -v DIScv/test3 -d formatted_data_GC.h5py -f rna   -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test1345 -v DIScv/test2 -d formatted_data_GC.h5py -f rna  -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer.py -t DIScv/test2345 -v DIScv/test1 -d formatted_data_GC.h5py -f rna   -bs 10

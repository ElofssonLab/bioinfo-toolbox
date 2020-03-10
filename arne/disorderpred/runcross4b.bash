#!/bin/bash -x

# Train Test Val
# 123 4 5
# 125 3 4
# 145 2 3
# 345 1 2
# 234 5 1


/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test123.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test125.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test145.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro.h5py -f pro   -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test234.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro.h5py -f pro   -bs 10  



#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test123.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10  
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test125.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test145.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic   -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test234.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic   -bs 10  

/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test123.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro_GC.h5py -f pro -gc  -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test125.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gc  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test145.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro_GC.h5py -f pro -gc   -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro_GC.h5py -f pro -gc  -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test234.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro_GC.h5py -f pro -gc   -bs 10  

/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test123.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gcgenomic -gc -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test125.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gcgenomic -gc -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test145.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gcgenomic -gc  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gcgenomic -gc -bs 10
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test234.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gcgenomic -gc  -bs 10  




#!/bin/bash -x

#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10  
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gcgenomic  -bs 10  
#
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gc  -bs 10  
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gc  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gc  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gc  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro_GCgenomic.h5py -f pro -gc  -bs 10  
#
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10  
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro.h5py -f pro  -bs 10  
#
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gc -gcgenomic   -bs 10  
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gc -gcgenomic   -bs 10  
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gc -gcgenomic   -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gc -gcgenomic   -bs 10 
#/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_pro_GC_GCgenomic.h5py -f pro -gc -gcgenomic   -bs 10  

/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gcgenomic  -bs 10  &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gcgenomic  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gcgenomic  -bs 10 &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gcgenomic  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gcgenomic  -bs 10  &

/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gc  -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gc  -bs 10 &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gc  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gc  -bs 10 &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_rna_GCgenomic.h5py -f rna -gc  -bs 10  

/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_rna.h5py -f rna  -bs 10  &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_rna.h5py -f rna  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_rna.h5py -f rna  -bs 10 &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_rna.h5py -f rna  -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_rna.h5py -f rna  -bs 10  &

/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1234.sorted.GCgenomic -v DIScv/test5.sorted.GCgenomic -d formatted_data_rna_GC_GCgenomic.h5py -f rna -gc -gcgenomic   -bs 10  
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1235.sorted.GCgenomic -v DIScv/test4.sorted.GCgenomic -d formatted_data_rna_GC_GCgenomic.h5py -f rna -gc -gcgenomic   -bs 10  &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1245.sorted.GCgenomic -v DIScv/test3.sorted.GCgenomic -d formatted_data_rna_GC_GCgenomic.h5py -f rna -gc -gcgenomic   -bs 10 
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test1345.sorted.GCgenomic -v DIScv/test2.sorted.GCgenomic -d formatted_data_rna_GC_GCgenomic.h5py -f rna -gc -gcgenomic   -bs 10 &
/opt/singularity3/bin/singularity exec --nv tf_image.sif python3 bin/dis_trainer3.py -t DIScv/test2345.sorted.GCgenomic -v DIScv/test1.sorted.GCgenomic -d formatted_data_rna_GC_GCgenomic.h5py -f rna -gc -gcgenomic   -bs 10  


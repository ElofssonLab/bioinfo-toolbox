#!/bin/bash -x
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f pro -gc -kingdom &
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f rna -gc -kingdom 
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f pro -gc 
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f rna -gc 
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f pro  
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f rna

/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -gc -gcgenomic -f pro -d ./
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -gcgenomic -f pro -d ./
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -gc -f pro -d ./
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -f pro -d ./

#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -final -gc -gcgenomic -f pro -d ./
#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -final -gcgenomic -f pro -d ./
#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -final -gc -f pro -d ./
#/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross2.py -final -f pro -d ./


/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross4.py  -gc -gcgenomic -f pro -d ./
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross4.py  -gcgenomic -f pro -d ./
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross4.py  -gc -f pro -d ./
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross4.py  -f pro -d ./

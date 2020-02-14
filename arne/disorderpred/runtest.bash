#!/bin/bash -x
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f pro -gc -kingdom &
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f rna -gc -kingdom &
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f pro -gc 
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f rna -gc &
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f pro  &
/opt/singularity3/bin/singularity exec --nv tf_sandbox/ python3 bin/dis_test_cross.py  -d ./ -f rna

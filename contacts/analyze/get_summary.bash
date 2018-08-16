#!/bin/bash -x


# remi
for i in `grep -cE "fa_.*fa_" PF*/*csv | grep -v :0 | sed s/:1//g ` ; do rm $i ; done


echo "target ,  model ,  ali ,  num ,  mindist ,  maxdist ,  length ,  TM ,  Pcons ,  cns ,  noe ,  ProQ2D ,  ProQ3D , Pcons-all, Meff " > summary-all.csv
cat PF*/*cm_summary.csv | grep -v target >> summary-all.csv

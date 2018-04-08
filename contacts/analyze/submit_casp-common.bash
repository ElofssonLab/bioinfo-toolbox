#!/bin/bash -x

# This is a very simple script that just submits one PconsC4 prediction and the hightest ranked COnFOLD model


MSA=JH0.001
# Submit RR

mail models@predictioncenter.org -s "Submission from Elofsson" -r arne@bioinfo.se < $1/$1.seq.${MSA}.rr

# Fix and submit TS



bin2/formatpdb.bash $1/$1.seq.${MSA}.rr_2.5_cm/stage1/$1.seq_1.pdb > $1/$1.pdb

mail models@predictioncenter.org -s "Submission from Elofsson" -r arne@bioinfo.se < $1/$1.pdb

#!/bin/bash


echo -n "PFRMAT TS
TARGET "
grep FILENAME $1 | sed 's/.*=\"//g' | sed 's/.seq.*//g'
echo "AUTHOR 5229-7541-3942
REMARK GROUP Elofsson
REMARK DEVELOPER Mirco Michel, Davide Menendez-Hurtado and Arne Elofsson 
METHOD PconsC4 - PconsFold
MODEL  1
PARENT N/A"
grep -v "END" $1
echo "TER"
echo "END"


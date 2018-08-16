#!/bin/bash



echo -n "PFRMAT RR
TARGET "
grep \> $1 | sed s/\>//g
echo "AUTHOR 5229-7541-3942
REMARK GROUP Elofsson
REMARK DEVELOPER Mirco Michel, Davide Menendez-Hurtado and Arne Elofsson 
METHOD PconsC4
MODEL  1"
grep -v \> $1
head -30000 $2 | gawk '{print $1, $2, 0, 8, $3 }' 
echo "END"


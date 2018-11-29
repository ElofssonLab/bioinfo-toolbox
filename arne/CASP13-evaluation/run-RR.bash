#!/bin/bash

grep --no-filename RR *L.L2.txt | sed s/.*RR/RR/g | sed s/-D[0-9]//g | sort > allresults.txt

for i in `gawk '{print $1}' hardtargets.txt`
do
    grep --no-filename RR $i*.L.L2.txt | sed s/.*RR/RR/g | sed s/-D[0-9]//g
done  | sort > hardresults.txt


gawk '{print $1}' allresults.txt | sort -u > methods.txt
gawk '{print $1}' hardresults.txt | sort -u > hardmethods.txt

for i in `cat methods.txt`
do
    echo -n $i " "
    grep $i allresults.txt | gawk '{i++; sum+=$3};END{print i,sum/i}'
done | sort -rnk 3 | sed s/RR152_1/PconsC4/g | sed s/RR497_1/GaussDCA/g | sed s/RR498_1/RaptorX-Contact/g | sed s/RR032_1/TripletRes1/g  | sed s/RR180_1/TripletRes2/g  | sed s/RR323_1/TripletRes3/g | sed s/RR491_1/DMP/g| sed s/RR125_1/TripletRes4/g | sed s/RR189_1/TripletRes5/g | sed s/RR352_1/RRMD/g | sed s/RR106_1/RRMD2/g | sed s/RR224_1/Defini/g | sed s/RR036_1/Zhang-contact/g  > summary-f1.txt

for i in `cat hardmethods.txt`
do
    echo -n $i " "
    grep $i hardresults.txt | gawk '{i++; sum+=$4};END{print i,sum/i}'
done | sort -rnk 3 | sed s/RR152_1/PconsC4/g | sed s/RR497_1/GaussDCA/g | sed s/RR498_1/RaptorX-Contact/g | sed s/RR032_1/TripletRes1/g  | sed s/RR180_1/TripletRes2/g  | sed s/RR323_1/TripletRes3/g | sed s/RR491_1/DMP/g | sed s/RR125_1/TripletRes4/g | sed s/RR189_1/TripletRes5/g | sed s/RR352_1/RRMD/g | sed s/RR106_1/RRMD2/g | sed s/RR224_1/Defini/g | sed s/RR036_1/Zhang-contact/g > summary-f1-hard.txt


for i in `cat methods.txt`
do
    echo -n $i " "
    grep $i allresults.txt | gawk '{i++; sum+=$4};END{print i,sum/i}'
done | sort -rnk 3 | sed s/RR152_1/PconsC4/g | sed s/RR497_1/GaussDCA/g | sed s/RR498_1/RaptorX-Contact/g | sed s/RR032_1/TripletRes1/g  | sed s/RR180_1/TripletRes2/g  | sed s/RR323_1/TripletRes3/g | sed s/RR491_1/DMP/g| sed s/RR125_1/TripletRes4/g | sed s/RR189_1/TripletRes5/g | sed s/RR352_1/RRMD/g | sed s/RR106_1/RRMD2/g | sed s/RR224_1/Defini/g  | sed s/RR036_1/Zhang-contact/g > summary-ppv.txt

for i in `cat hardmethods.txt`
do
    echo -n $i " "
    grep $i hardresults.txt | gawk '{i++; sum+=$4};END{print i,sum/i}'
done | sort -rnk 3 | sed s/RR152_1/PconsC4/g | sed s/RR497_1/GaussDCA/g | sed s/RR498_1/RaptorX-Contact/g | sed s/RR032_1/TripletRes1/g  | sed s/RR180_1/TripletRes2/g  | sed s/RR323_1/TripletRes3/g | sed s/RR491_1/DMP/g | sed s/RR125_1/TripletRes4/g | sed s/RR189_1/TripletRes5/g | sed s/RR352_1/RRMD/g | sed s/RR106_1/RRMD2/g | sed s/RR224_1/Defini/g | sed s/RR036_1/Zhang-contact/g > summary-ppv-hard.txt



# sed s/RR498_1/RaptorX-Contact/g | sed s/RR032_1/TripletRes1/g  | sed s/RR180_1/TripletRes2/g /g | sed s/RR323_1/TripletRes3/g | sed s/RR491_1/DMP/g| sed s/RR125_1/TripletRes4/g /g| sed s/RR189_1/TripletRes5/g | sed s/RR352_1/RRMD/g | sed s/RR106_1/RRMD2/g | sed s/RR224_1/Defini/g

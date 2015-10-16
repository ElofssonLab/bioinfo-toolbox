#!/bin/sh

#2009-07-15 
# runpsiblast

hosts="
grisman-01
grisman-02
grisman-03
grisman-04
grisman-05
grisman-06
grisman-07
grisman-08
grisman-09
grisman-10
grisman-11
grisman-12
grisman-13
grisman-14
grisman-15
grisman-16
grisman-17
grisman-18
grisman-19
grisman-20
grisman-21
grisman-22
grisman-23
grisman-24
grisman-25
grisman-26
grisman-27
grisman-28
grisman-29
grisman-30
grisman-31
grisman-32
grisman-33
grisman-34
grisman-35
linux12
linux13
linux14
linux15
linux16
linux17
linux18
linux19
linux20
linux21
grisman-17
grisman-18
grisman-19
grisman-20
grisman-21
grisman-22
grisman-23
grisman-24
grisman-25
grisman-26
grisman-27
grisman-28
grisman-29
grisman-30
grisman-31
grisman-32
grisman-33
grisman-34
grisman-35
linux12
linux13
linux14
linux15
linux16
linux17
linux18
linux19
linux20
linux21
grisman-17
grisman-18
grisman-19
grisman-20
grisman-21
grisman-22
grisman-23
grisman-24
grisman-25
grisman-26
grisman-27
grisman-28
grisman-29
grisman-30
grisman-31
grisman-32
grisman-33
grisman-34
grisman-35
linux12
linux13
linux14
linux15
linux16
linux17
linux18
linux19
linux20
linux21
grisman-17
grisman-18
grisman-19
grisman-20
grisman-21
grisman-22
grisman-23
grisman-24
grisman-25
grisman-26
grisman-27
grisman-28
grisman-29
grisman-30
grisman-31
grisman-32
grisman-33
grisman-34
grisman-35
linux12
linux13
linux14
linux15
linux16
linux17
linux18
linux19
linux20
linux21
grisman-17
grisman-18
grisman-19
grisman-20
grisman-21
grisman-22
grisman-23
grisman-24
grisman-25
grisman-26
grisman-27
grisman-28
grisman-29
grisman-30
grisman-31
grisman-32
grisman-33
grisman-34
grisman-35
linux12
linux13
linux14
linux15
linux16
linux17
linux18
linux19
linux20
linux21
"

#for test set
pdbaapath=nrpdb99aa
idlistFile=newTest.variedEvalue.idlist
outpath=newTestVariedEvalue

mkdir -p $outpath

N=2421
Nbeg=0
Nend=2440
Nstep=50

((cntSplitID=0))
for ((i=Nbeg; i <=Nend; i +=Nstep)); do
    start[$cntSplitID]=$i
    end[$cntSplitID]=`expr $i + $Nstep`
    ((cntSplitID++))
done
numSplitID=$cntSplitID
echo "numSplitID = $numSplitID"

((cntSplitID=0))
for host in ${hosts}; do
    ssh nanjiang@$host <<EOF 
    screen -r
    cd  /usr/project/xtmp/tmp/nj/wk/psipredTestNew
    echo $host
    echo "./runpsiblast-variedEvalue.sh --pdbaa $pdbaapath -l $idlistFile -a 2 -b ${start[$cntSplitID]} -e ${end[$cntSplitID]} --outpath $outpath > /dev/null 2>&1 &"
    ./runpsiblast-variedEvalue.sh --pdbaa $pdbaapath -l $idlistFile -a 2 -b ${start[$cntSplitID]} -e ${end[$cntSplitID]} --outpath $outpath > /dev/null 2>&1 &
    screen -X detach
    exit
EOF
   break
    ((cntSplitID++))
    if [ $cntSplitID -ge $numSplitID ] ; then
        break
    fi
done

#!/bin/bash
# runPsiblast_para.sh
# 2010-08-11

usage="
Usage:     runPsiblast_para.sh argument-list
Options:
    -l|--list file          : idlistFile
    -aapath|--aapath path   : set the path for sequence files for ids
    -outpath|--outpath path : set the output path
    -workdir|--workdir path : set WORKINGDIR
    -sleep <int>            : set sleep time in seconds, this is only for submitting batch jobs on multiple hosts
                            : to wait a certain amount of time, so that the login will not be blocked
                            : when the program starts
    -j int                  : num of iterations for blastpgp, default=4
    -db       file          : set the blastdb, default=
    -bb       int           : __blastpgp_num_align_b
    -h|--help               : print this help message and exit

Created 2009-12-15, updated 2010-08-11, Nanjiang
"
function PrintHelp()
{
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

hostSet=( #{{{
grisman-17   5
grisman-18   5
grisman-19   5
grisman-21   5
grisman-22   5
grisman-23   5
grisman-25   5
grisman-26   5
grisman-27   5
grisman-28   5
grisman-31   5
grisman-32   5
grisman-33   5
grisman-34   5
linux22      12
linux23      12
linux24      12
linux25      12
linux26      12
linux28      12
linux29      12
linux30      12
)      #}}}


numHost=${#hostSet[*]}
((numHost/=2))
echo "numHost=$numHost"

sumCPUs=0
for ((i = 0; i< numHost; i ++)); do
    ((m=i*2))
    ((n=m+1))
    hosts[$i]=${hostSet[$m]}
    CPUs[$i]=${hostSet[$n]}
    ((sumCPUs+=CPUs[$i]))
    #echo "${hosts[$i]}  ${CPUs[$i]}"
done

echo $sumCPUs

idListFile=
aaPath=
WORKINGDIR=
outpath=

runPsiblastOption=
isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -l|--list) idListFile=$2;shift;;
            -aapath|--aapath) aaPath=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -workdir|--workdir) WORKINGDIR=$2;shift;;
            -sleep|--sleep) runPsiblastOption="$runPsiblastOption -sleep $2"; shift;;
            -j) runPsiblastOption="$runPsiblastOption -j $2"; shift;;
            -db) runPsiblastOption="$runPsiblastOption -db $2"; shift;;
            -bb) runPsiblastOption="$runPsiblastOption -bb $2"; shift;;
            -*) echo "Error! Wrong argument: $1"; exit;;
        esac
    else
        echo "Error! Wrong argument: $1";
        exit
    fi
    shift
done

if [ ! -f "$idListFile" ]; then 
    echo "Error! idListFile = \"$idListFile\" does not exist. Exit..."
    exit
fi

if [ ! -d "$aaPath" ]; then 
    echo "Error! aaPath = \"$aaPath\" does not exist. Exit..."
    exit
fi
if [ "$WORKINGDIR"  == "" ]; then 
    echo "Error! WORKINGDIR = \"$WORKINGDIR\" does not set. Exit..."
    exit
fi

if [ "$outpath" == "" ]; then 
    outpath=$WORKINGDIR
fi

if [ "${WORKINGDIR:0:1}" != "/" ];then
    WORKINGDIR=$PWD/$WORKINGDIR # add the absolut path
fi
if [ "${idListFile:0:1}" != "/" ];then
    idListFile=$PWD/$idListFile # add the absolut path
fi
if [ "${aaPath:0:1}" != "/" ];then
    aaPath=$PWD/$aaPath # add the absolut path
fi
if [ "${outpath:0:1}" != "/" ];then
    outpath=$PWD/$outpath # add the absolut path
fi


mkdir -p $outpath
mkdir -p $WORKINGDIR

# Get the number of IDs to run on each host
# numID4Host[i] = (N / sumCPUs) * CPUs[i]

N=`wc -l $idListFile | awk '{print $1}' `
Ratio_N_SUMCPU=`e "$N / $sumCPUs "` 
# | awk '{printf "%.0f", $1}'` 
for ((i=0;i<numHost;i++));do
    numID4Host[$i]=` e  "$Ratio_N_SUMCPU * ${CPUs[$i]} + 1 " | awk '{printf "%.0f", $1}'`
#    echo ${numID4Host[$i]}
done


Nbeg=0
Nend=$N
((cntSplitID=0))
((accCPUs=0))
for ((i=0; i <numHost; i ++)); do
      begin[$cntSplitID]=$accCPUs
      ((end[$cntSplitID] = accCPUs+numID4Host[$i]))
      ((accCPUs += numID4Host[$i]))
      ((cntSplitID++))
#      echo "i=$i, begin = ${begin[$i]}, end = ${end[$i]}, numID4Host = ${numID4Host[$i]}"
      if [ $accCPUs -ge $Nend ]; then 
          break
      fi
done
numSplitID=$cntSplitID
echo "numSplitID = $numSplitID"

for ((cntSplitID=0;cntSplitID<numSplitID;cntSplitID++));do
    BEG=${begin[${cntSplitID}]}
    END=${end[${cntSplitID}]}
    HOST=${hosts[${cntSplitID}]}
    errLogFile=log/run_predzinc_para.$HOST.$$.err
    #============  Remote Code =====================
    ssh nanjiang@$HOST <<EOF 
    #!/bin/bash
    echo $HOST
    screen -r
    cd $WORKINGDIR
    mkdir -p log
    echo "runpsiblast-variedEvalue.sh $runPsiblastOption -l $idListFile --pdbaa $aaPath -b $BEG -e $END -a ${CPUs[$cntSplitID]} -outpath $outpath &> $errLogFile &"
    runpsiblast-variedEvalue.sh $runPsiblastOption -l $idListFile --pdbaa $aaPath -b $BEG -e $END -a ${CPUs[$cntSplitID]} -outpath $outpath &> $errLogFile &

    echo "errLogFile = $errLogFile"
    echo
    screen -X detach
    exit
EOF
    #============  Remote Code =====================
done


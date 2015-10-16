#!/bin/bash
# submit jobs on uppmax
# Filename:   submit_run_hhsearch_pairwise_uppmax.sh
nodename=`uname -n`


PrintHelp() {
    echo "$usage"
}
binpath=$DATADIR3/wk/MPTopo/src

case $nodename in 
    *uppmax.uu.se) echo "OK, you are running from uppmax node" ;;
    *) "Error. This script should be run from uppmax login node. Exit"; exit 1;;
esac

progname=`basename $0`
usage="
usage: $progname tableinfofile [-outpath outpath]

-project STR    Set project, (default: snic2014-8-12)
-hhprofile STR  Set path to hhprofile
-num     INT    Number per split, (default: 2000)
-time    STR    Specify time to run, default: 24:00:00
                e.g. 
                1-20:00:00
                2-20:00:00
"
infile=
outpath=./
hhprofile=
topofile=
project=
numseqperjob=2000
timetorun=24:00:00 


isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        infile=$1
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) outpath=$2;shift;;
            -project|--project) project=$2;shift;;
            -num|--num) numseqperjob=$2;shift;;
            -topo|--topo) topofile=$2;shift;;
            -time|--time) timetorun=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

if [ "$hhprofile" == "" ];then
    echo "hhprofile not set. exit"
    exit 1
fi
if [ "$topofile" == "" ];then
    echo "topofile not set. exit"
    exit 1
fi

if [ "$infile" == ""  ]; then 
    echo "infile not set. Exit" >&2
    exit 1
elif [ ! -s "$infile" ] ; then
    echo "infile '$infile' not exist or empty. Exit. " >&2
else 
    echo "Command line: $commandline"
    echo "Start at    `/bin/date`"
    res1=$(/bin/date +%s.%N)
    basename=`basename $infile`
    rootname=${basename%.*}

    numpair=`cat  $infile | wc -l `
    if [ "$numpair" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        exit 1
    elif [ $numpair -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        exit 1
    else
        mkdir -p $outpath
        mkdir -p $outpath/log
        mkdir -p $outpath/script
        #seqid2resultpathMapFile=$outpath/id2pathmap.txt
        #cat /dev/null > $seqid2resultpathMapFile
        echo "$numpair pairs are going to be processed by hhalign_pairwise"
        splittedpath=$outpath/splitted
        splitline.sh $infile -mline $numseqperjob -outpath $splittedpath
        numsplitted=`/bin/find $splittedpath -name "*split*" | wc -l `
        for ((i=0;i<=numsplitted;i++)); do
            splittedpairlistfile=$splittedpath/${rootname}.split${i}
            logfile=$outpath/log/${rootname}_split${i}.log
            jobscriptfile=$outpath/script/run_hhalign_pairwise_split${i}.sh
            jobname=hhalign_sp${i}
            suboutpath=$outpath/split_$i
            mkdir -p $suboutpath
            if [ -s "$splittedpairlistfile" ]; then 
                echo "\
#!/bin/bash
$DATADIR3/wk/MPTopo/src/run_hhsearch_pairwise.py -l $splittedpairlistfile -hhprofile $hhprofile -topofile $topofile -outpath $suboutpath 2> $logfile 1> /dev/null" > $jobscriptfile
                echo sbatch -A $project  -p core -n 1 -t $timetorun -J $jobname $jobscriptfile
                sbatch -A $project  -p core -n 1 -t $timetorun -J $jobname $jobscriptfile
                #python $DATADIR3/wk/MPTopo/bin/getfastaid.py $splittedfastafile | awk -v dir=split_$i '{print $1 "\t" dir}' >> $seqid2resultpathMapFile
            fi
        done
        echo "Finished at `/bin/date`"
    fi
fi

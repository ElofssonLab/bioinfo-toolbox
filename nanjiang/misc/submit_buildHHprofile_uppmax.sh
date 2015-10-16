#!/bin/bash
# submit jobs on uppmax
# it takes about 2 minutes to run one sequence with one core
nodename=`uname -n`
rundir=`dirname $0`

binpath=$rundir

case $nodename in 
    *uppmax.uu.se) echo "OK, you are running from uppmax node" ;;
    *) "Error. This script should be run from uppmax login node. Exit"; exit 1;;
esac

usage="
usage: submit_buildHHprofile_uppmax.sh fasta-seq-file [-outpath outpath]

-project STR    Set project, (default: snic2013-11-2)
-num     INT    Number per split, (default: 200)
-time    STR    Specify time to run, default: 24:00:00
                e.g. 
                1-20:00:00
                2-20:00:00
"
infile=
outpath=./
hhdb=$HHDB/uniprot20_02Sep11
project=snic2013-11-2
numseqperjob=200
timetorun=24:00:00 


PrintHelp(){
    echo "$usage"
}

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
            -hhdb|--hhdb) hhdb=$2;shift;;
            -num|--num) numseqperjob=$2;shift;;
            -time|--time) timetorun=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1">&2; exit;;
        esac
    else
        infile=$1
    fi
    shift
done

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

    numseq=`python $binpath/countseq.py $infile -nf `
    if [ "$numseq" == "" ]; then 
        echo "Fatal. countseq.py error. Exit."  >&2
        exit 1
    elif [ $numseq -le 0 ] ; then
        echo "Zero sequence found in file '$infile'. Exit."  >&2
        exit 1
    else
        mkdir -p $outpath
        mkdir -p $outpath/log
        mkdir -p $outpath/script
        seqid2resultpathMapFile=$outpath/id2pathmap.txt
        cat /dev/null > $seqid2resultpathMapFile
        echo "$numseq sequences are going to be processed by buildHHprofile"
        #numseqperjob=` expr  200 `
        # about 3 minutes per protein
        splittedpath=$outpath/splitted
        python $binpath/splitfasta.py -i $infile -nseq $numseqperjob -outpath $splittedpath -ext fa
        basename_infile=`basename $infile`
        rootname_infile=${basename_infile%.*}
        numsplitted=`/bin/find $splittedpath -name "${rootname_infile}*.fa" | wc -l `
        for ((i=0;i<=numsplitted;i++)); do
            splittedfastafile=$splittedpath/${rootname}_${i}.fa
            logfile=$outpath/log/${rootname}_split${i}.log
            jobscriptfile=$outpath/script/run_buildHHprofile_split${i}.sh
            jobname=hhblits_sp${i}
            suboutpath=$outpath/split_$i
            mkdir -p $suboutpath
            if [ -s "$splittedfastafile" ]; then 
                numfinished=`find $suboutpath -name "*.hhm" | wc -l`
                subnumseq=`python $binpath/countseq.py $splittedfastafile -nf `
                if [ "$numfinished" -eq "$subnumseq" ]; then
                    continue
                fi
                echo "\
#!/bin/bash
$DATADIR3/wk/MPTopo/bin/buildHHprofile.sh $splittedfastafile -no-overwrite -outpath $suboutpath -hhdb $hhdb > $logfile 2>&1" > $jobscriptfile
                echo sbatch -A $project  -p core -n 1 -t $timetorun -J $jobname $jobscriptfile
                sbatch -A "$project"  -p core -n 1 -t $timetorun -J $jobname $jobscriptfile
                python $DATADIR3/wk/MPTopo/bin/getfastaid.py $splittedfastafile | awk -v dir=split_$i '{print $1 "\t" dir}' >> $seqid2resultpathMapFile
            fi
        done
        echo "Finished at `/bin/date`"
    fi
fi

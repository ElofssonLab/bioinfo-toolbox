#!/bin/bash

# kalign
# muscle
# clustalw
# clustalo
# t_coffee
# probcons
# run_kalignP.sh
# mafft
# probcons

# Note: t_coffee runs on multiple cores, the real time used is not comparable to
# others

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 
usage="
Usage:  $progname [-l LISTFILE] [-outpath DIR]  [-q]
        $wspace [FILE [FILE ...]] 
        $wspace [-prog PROG1 -prog PROG2]

Description: Given a file of multiple sequences and run multiple sequence
             alignments for all methods in the program list
Options:
  -outpath DIR      Set outpath, (default: the same as input file)
  -l       FILE     Set the fileListFile, one filename per line
  -prog     STR     Set the program name in command line.
  -proglist  FILE   Set the program list file, one item per line
  -showtime y|n     Whether show time, (default: yes)
  -overwrite y|n    Whether overwrite the result, (default: no)
  -x, -compress y|n Whether compress the msa file, (default: no)
  -h, --help        Print this help message and exit

Created 2013-02-13, updated 2013-02-21, Nanjiang Shu 

Examples:
#run multiple msa and compress the result file
  $program -l test.fa -x yes

"

nodename=`uname -n`
TMPDIR=/scratch/$$
case $nodename in 
    *pdc.kth.se) 
        datadir3=$DATADIR3 
        TMPDIR=/scratch/$$
        ;;
    *uppmax.uu.se|triolith*)
        datadir3=$DATADIR3
        if [ $SNIC_TMP ]; then
            TMPDIR=$SNIC_TMP
        else
            TMPDIR=/tmp/$$
        fi
        ;;
    *vault*) 
        datadir3=/big/nanjiang/data3
        TMPDIR=$datadir3/tmp/runmsa_$$
        ;;
    *) 
        datadir3=/data3 
        if [ -d /scratch ]; then 
            TMPDIR=/scratch/$$
        else
            TMPDIR=/tmp/$$
        fi
        ;;
esac
if [ ! -d $TMPDIR ] ; then 
    mkdir -p $TMPDIR
fi

path_kalignP=$datadir3/program/run_kalignP
path_mafft=$datadir3/usr/bin


PrintHelp(){
    echo "$usage"
}
IsProgExist(){ #{{{
    # usage: IsProgExist prog
    # prog can be both with or without absolute path
    type -P $1 &>/dev/null \
        || { echo The program \'$1\' is required but not installed. \
        Aborting $0 >&2; isAllProgramExist=0; return 1; }
    return 0
}
#}}}
in_array() {
    local hay needle=$1
    shift
    for hay; do
        if [ "$hay" == "$needle" ]; then
            echo 1
            return 0
        fi
    done
    echo 0
    return 1
}

CheckMSAProgram(){ #{{{
    for((i=0;i<numProg;i++)); do
        local prog=${progList[$i]}
        case $prog in
            kalignP)
                IsProgExist $path_kalignP/run_kalignP.sh 
                ;;
            kalign) 
                IsProgExist kalign
                ;;
            muscle*) 
                IsProgExist muscle
                ;;
            clustalw)
                IsProgExist clustalw
                ;;
            clustalo*)
                IsProgExist clustalo
                ;;
            t_coffee)
                IsProgExist t_coffee
                ;;
            probcons)
                IsProgExist probcons
                ;;
            mafft*)
                IsProgExist $path_mafft/mafft
                ;;
        esac
    done
}
#}}}

CompressMSAFile(){ #$seqfile $msafile#{{{
    local seqfile=$1
    local msafile=$2
    if [ ! -s $msafile ]; then
        /bin/rm -f $msafile
        echo Delete empty mfa file $msafile
        return
    else
        num_aaseq=`grep "^>" $seqfile | wc -l`
        num_msaseq=`grep "^>" $msafile | wc -l`
        numdifferentseqlength=`python $exe_getseqlen $msafile | sort -u | wc -l`

        if [ "$num_aaseq" == "$num_msaseq" -a "$numdifferentseqlength" == "1" ]; then
            echo "gzip -fN $msafile"
            gzip -fN $msafile
        else
            echo MSA file $msafile incomplete, Delete $msafile.
            /bin/rm -f $msafile
        fi
    fi
}
#}}}
RunMSA(){ # seqfile prog $outpath#{{{
    # output format: fasta
    local res1=
    local res2=
    if [ $isShowRunningTime -eq 1 ]; then
        res1=$(/bin/date +%s.%N)
    fi
    local seqfile=$1
    local prog=$2
    local outpath=$3
    local numseq=`grep "^>" $seqfile | wc -l`

    local basename=`basename "$seqfile"`
    local rootname=${basename%.*}

    local outfile=$outpath/$rootname.$prog.mfa

    (case $prog in
        kalignP)
            $path_kalignP/run_kalignP.sh $seqfile -o $outfile \
                -f fasta -q -tmp $TMPDIR
            ;;
        kalign) 
            kalign $seqfile -q -o $outfile -f fasta
            ;;
        muscle) 
            # default: fasta format
            muscle -in $seqfile -quiet -out $outfile
            ;;
        muscle_maxiter2) 
            # default: fasta format, for higher speed, set maximum iteration
            # times to 2
            muscle -in $seqfile -quiet -out $outfile -maxiters 2
            ;;
        clustalw)
            clustalw $seqfile -outfile=$outfile -output=FASTA -quiet >/dev/null
            ;;
        clustalo)
            clustalo -i $seqfile -o $outfile --outfmt FASTA --iter=2 --force --threads=1
            ;;
        clustalo_auto)
            clustalo -i $seqfile -o $outfile --outfmt FASTA \
                --auto --force --threads=1
            ;;
        clustalo*)
            iter_clustalo=`echo $prog | sed 's/clustalo//g'`
            clustalo -i $seqfile -o $outfile --outfmt FASTA --iter=$iter_clustalo --force --threads=1
            ;;
        t_coffee)
            t_coffee $seqfile -outfile $outfile -output fasta -quiet
            ;;
        probcons)
            probcons $seqfile > $outfile
            ;;
        mafft)
            $path_mafft/mafft --quiet $seqfile > $outfile
            ;;
        mafft_auto)
            $path_mafft/mafft --auto --quiet $seqfile > $outfile
            ;;
        *)
            return 1
            ;;
    esac
    ) 2> /dev/null
    if [ $isShowRunningTime -eq 1 ]; then
        res2=$(/bin/date +%s.%N)
    fi
    printf "Alignment of %s with %d sequences by %s costs %.6F\n" \
        "$seqfile" $numseq "$prog" $(echo "$res2 - $res1"|/usr/bin/bc)
}
#}}}
Exec_RunMSA(){ #seqfile prog #{{{

    local seqfile="$1"
    local prog="$2"

    if [ ! -s "$seqfile" ]; then
        echo "seqfile $seqfile does not exist or empty. Ignore"
        return 1
    fi

    local outpath=
    if [ "$g_outpath" != "" ];then
        outpath=$g_outpath
    else
        outpath=`dirname $seqfile`
    fi

    local num_aaseq=`grep "^>" $seqfile | wc -l`
    local basename=`basename "$seqfile"`
    local rootname=${basename%.*}

    local msafile=$outpath/$rootname.$prog.mfa
    local gzfile=$outpath/$rootname.$prog.mfa.gz
    local isResultExist=0

    if [ -s $gzfile ]; then
        isResultExist=1
    elif [ -s $msafile ];then
        local num_msaseq=`grep "^>" $msafile | wc -l`
        local numdifferentseqlength=`python $exe_getseqlen $msafile | sort -u | wc -l`
        if [ "$num_aaseq" == "$num_msaseq" -a "$numdifferentseqlength" == "1" ]; then
            isResultExist=1
        fi
    fi

    if [ $isResultExist -eq 1 -a $isOverWrite -eq 0 ]; then
        echo Result $msafile already exist. Ignore
        return 
    fi
    local msg=$((time -p RunMSA "$seqfile" "$prog" "$outpath" ) 2>&1 )
    echo $msg

    if [ $isCompressMSAFile -eq 1 ];then
        CompressMSAFile $seqfile $msafile
    fi
}
#}}}

isQuiet=0
g_outpath=
fileListFile=
fileList=()
isShowRunningTime=1
isOverWrite=0
progListFile=

default_progList=(
"kalign"
"kalignP"
"muscle_maxiter2"
"mafft"
"mafft_auto"
"clustalo"
"clustalo_auto"
"clustalw"
"muscle"
"t_coffee"
"probcons"
)

optProgList=()


if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isCompressMSAFile=0
isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        fileList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -outpath|--outpath) g_outpath=$2;shift;;
            -o|--o) outfile=$2;shift;;
            -proglist|--proglist) progListFile=$2;shift;;
            -prog|--prog) optProgList+=("$2");shift;;
            -l|--l|-list|--list) fileListFile=$2;shift;;
            -st|--st|-showtime|--showtime)
                opt=$2
                if [ "${2:0:1}" == "y" -o "${2:0:1}" == "Y" ]; then
                    isShowRunningTime=1
                else
                    isShowRunningTime=0
                fi
                shift;
                ;;
            -x|--x|-compress|--compress)
                opt=$2
                if [ "${2:0:1}" == "y" -o "${2:0:1}" == "Y" ]; then
                    isCompressMSAFile=1
                else
                    isCompressMSAFile=0
                fi
                shift;
                ;;
            -overwrite|--overwrite)
                opt=$2
                if [ "${2:0:1}" == "y" -o "${2:0:1}" == "Y" ]; then
                    isOverWrite=1
                else
                    isOverWrite=0
                fi
                shift;
                ;;
            -q|-quiet|--quiet) isQuiet=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else
        fileList+=("$1")
    fi
    shift
done

if [ "$progListFile" != ""  ]; then 
    if [ -s "$progListFile" ]; then 
        progList=()
        while read line
        do
            progList+=("$line")
#           if [ $(in_array "$line" ${default_progList[@]}) -eq 1 ]; then
#               progList+=("$line")
#           else
#               echo unrecognized program \"$line\". Ingore. >&2
#           fi
        done < "$progListFile"
    else
        echo proglistfile \'$progListFile\' does not exist or empty. >&2
    fi
elif [ ${#optProgList[@]} -ge 1  ]; then
    num=${#optProgList[@]}
    for((i=0;i<num;i++));do
        prog=${optProgList[$i]}
        progList+=("$prog")
#       if [ $(in_array "$prog" ${default_progList[@]}) -eq 1 ]; then
#           progList+=("$prog")
#       else
#           echo unrecognized program \"$prog\". Ignore. >&2
#       fi
    done
else
    progList=(${default_progList[@]})
fi

exe_getseqlen=`which getseqlen.py`

if [ "$exe_getseqlen" == "" ]; then
    echo getseqlen.py does not exist. Exit >&2
    exit 1
fi


if [ "$fileListFile" != ""  ]; then 
    if [ -s "$fileListFile" ]; then 
        while read line
        do
            fileList+=("$line")
        done < $fileListFile
    else
        echo listfile \'$fileListFile\' does not exist or empty. >&2
    fi
fi

if [ "$g_outpath" != ""  -a ! -d "$g_outpath" ]; then
    mkdir -p "$g_outpath"
fi

numFile=${#fileList[@]}
if [ $numFile -eq 0  ]; then
    echo Input not set! Exit. >&2
    exit 1
fi

numProg=${#progList[@]}
# for((i=0;i<numProg;i++)); do
#     prog=${progList[$i]}
#     echo $prog
# done

# check progList
isAllProgramExist=1
CheckMSAProgram
if [ $isAllProgramExist -eq 0 ]; then
    echo Some MSA programs do not exist. Exit.
    exit 1
fi

for ((i=0;i<numProg;i++));do
    prog=${progList[$i]}
    for ((j=0;j<numFile;j++));do
        seqfile=${fileList[$j]}
        Exec_RunMSA "$seqfile" "$prog"
    done
done

/bin/rm -rf $TMPDIR

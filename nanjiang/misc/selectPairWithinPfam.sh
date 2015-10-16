#!/bin/bash
# select pair within pfam family by sequence identity


progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 

rundir=`dirname $0`
anamode=0  #default anamode = 0
binpath=$rundir/../bin
dgpredpath=$rundir/../../../program/dgpred_standalone
newscampiscriptpath=$rundir/../../../program/newscampiscript
run_kalignP_path=$rundir/../../../program/run_kalignP
scampi_dir=$rundir/../../../share/scampi
modhmm_bin=$rundir/../../../share/modhmm/bin
emboss_bin=$rundir/../../../usr/share/emboss/bin
cdhit_bin=$rundir/../../../bin
clustalo_bin=$rundir/../../../bin
kalign_bin=$rundir/../../../bin
my_needlebin=$rundir/../../../program/my_needle
isPrintVerboseInfo= # in bash, if [ $istrue ]. the empty variable means false and any non-empty variable means true 
blastbin=$BLASTBIN
#blastdb=$BLASTDB/uniprotWG100Filtered
blastdb=$BLASTDB/uniprotWG100
#blastdb=$BLASTDB/swissprot

numCPU=`grep process /proc/cpuinfo | wc -l`
numCPU=2

usage="
Usage: $progname [-b INT] [-e INT] [-v] [-outpath]
       $wspace [-anamode 0   -dbnameseq STR ]
       $wspace [-datapath DIR]
       $wspace -l LISTFILE  or pfamid [pfamid...]
Description:
    select pairs with high sequence identity with Pfam families

Options:
    -seqidt  FLOAT    Sequence identity threshold for CD-HIT, 0-100, (default: 80)
    -l FILE           Set the file containing ids (rootname) of sequences
    -b INT            Begin index for pfamid, default = 0
    -e INT            End index for pfamid, default = INF
    -outpath  DIR     Location to output the result, (default: ./)
    -anamode  INT     Set the analysis mode, 0, 1, 2 or 3 (default: $anamode)
                      0. given sequences and multiple alignments
    -dbnameseq STR    Set the database name for sequences
    -datapath  DIR    Set the datapath containing sequence files \$pfamid.fa
    -v, --verbose     Print verbose information
    -nc               Do not clean
    -h, --help        Print this help message and exit

Created 2012-03-22, updated 2012-03-22, Nanjiang Shu 

Examples:
#select pairs with sequence identity >80% (cd-hit seqidt)
    selectPairWithinPfam.sh -l sel64.pfamidlist -outpath result -dbnameseq seqdb -seqidt 80
nanjiang@sbc.su.se
"

PrintHelp(){ #{{{ 
    echo "$usage"
}
#}}}
IsProgExist(){ #{{{
    # usage: IsProgExist prog
    # prog can be both with or without absolute path
    type -P $1 &>/dev/null \
        || { echo The program \"$1\" is required but it\'s not installed. \
        Aborting. >&2; exit 1; }
}
#}}}
IsPathExist(){ #{{{
    # supply the effective path of the program 
    if ! test -d "$1"; then
        echo Directory \'$1\' does not exist. Aborting. >&2
        exit 1
    fi
}
#}}}
IsAllSeqIDUnique(){ #{{{
    local file=$1
    if [ ! -s  "$file" ]; then 
        echo File \'$1\' is empty >&2
        return 1
    else
        python $binpath/getseqlen.py  $file  --printid | awk '{print $1}' \
            | sort | uniq -d
        return 0
    fi
}
#}}}

CleanTempFile(){
    local pfamid=$1
    if [ $isClean -eq 1 ]; then
        rm -f $outpath/$pfamid.fa
        rm -f $outpath/$pfamid.nr${c_seqidt}.fa*
        rm -f $outpath/$pfamid.seqdb*
        rm -f $outpath/$pfamid.c${c_seqidt}.pairlist
    fi
}

SelectPairWithinPfam(){ # $pfamid#{{{ 
    # Created 2012-03-22, updated 2012-03-22, Nanjiang Shu 
    local pfamid=$1
    cp -f $datapath/$pfamid.fa $outpath/
    local seqfile=$outpath/$pfamid.fa

    if [ ! -s "$seqfile" ]; then 
        echo seqfile \'$seqfile\' does not exist or empty. Ignore. >&2
        CleanTempFile $pfamid
        return 1
    else
        local numseq=`python $binpath/countseq.py $seqfile -nf `
        if [ $numseq -le 1 ]; then 
            echo Too few sequences \($numseq\) for pfamid $pfamid. At least 2 \
                sequences are needed. Ignore. >&2
            CleanTempFile $pfamid
            return 1
        fi
    fi

    ### Step 1: CD-hit to find clusters with high sequence identity
    ### using scampi using the scampi without xslt
    local nrseqfile=$outpath/$pfamid.nr${c_seqidt}.fa
    c_seqidt_float=`e $c_seqidt / 100`
    $cdhit_bin/cd-hit -i $seqfile -o $nrseqfile -c $c_seqidt_float -n 5 -M 2000
    clstrfile=$outpath/$pfamid.nr${c_seqidt}.fa.clstr 

    ### Step 2: Get pairs with seqidt above threshold
    if [ ! -s "$clstrfile" ]; then 
        echo clstrfile \'$clstrfile\' does not exist or empty. Ignore. >&2
        CleanTempFile $pfamid
        return 1
    else
        selectedpairfile=$outpath/$pfamid.c${c_seqidt}.pairlist
        python $binpath/cluster2pair_cdhit.py $clstrfile -o $selectedpairfile
    fi

    ### Step 3: align pairs by my_needle 
    if [ ! -s "$selectedpairfile" ]; then
        echo selectedpairfile \'$selectedpairfile\' does not exist or empty. Ignore. >&2
        CleanTempFile $pfamid
        return 1
    else
        outTableFile=$outpath/$pfamid.tableinfo.txt
        outAlnFile=$outpath/$pfamid.pairaln.fa
        tmp_seqdb=$outpath/$pfamid.seqdb
        $binpath/indexfasta.py $seqfile -dbname $tmp_seqdb
        $my_needlebin/my_needle -list $selectedpairfile -mode 1 -seqdb $tmp_seqdb -m 1 -table $outTableFile -o $outAlnFile > /dev/null
    fi
    CleanTempFile $pfamid
}
#}}}


if [ $# -lt 1 ]; then
    PrintHelp
    exit 1
fi

commandline="$0 $*"

errFile=$(mktemp /tmp/tmperr.$progname.XXXXXXX) \
    || { echo Failed to create temp file >&2; exit 1; }

idListFile=
dbnameseq=$rundir/../pfamAna/pfamfullseq.selTM_uniq
datapath=
outpath=
begin=0
end=9999999
isClean=1

idList=
c_seqidt=80


isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then 
        isNonOptionArg=0
        idList="$idList $1"
    elif [ "$1" == "--" ]; then
        isNonOptionArg=1
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -l|--l|--list|-list|--idlist) idListFile=$2;shift;;
            -b) begin=$2;shift ;;
            -e) end=$2;shift ;;
            -datapath|--datapath) datapath=$2;shift;;
            -dbnameseq|--dbnameseq) dbnameseq=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -seqidt|--seqidt) c_seqidt=$2;shift;;
            -anamode|--anamode) anamode=$2;shift;;
            -nc|--nc) isClean=0;;
            -v|-verbose|--verbose) isPrintVerboseInfo=1;shift;;
            -*) echo Error! Wrong argument: $1; exit;;
        esac
    else
        idList="$idList $1"
    fi
    shift
done

IsProgExist /usr/bin/bc
IsProgExist /bin/date
IsProgExist $binpath/my_extractdb.py
IsProgExist $binpath/countseq.py
IsProgExist $binpath/e

if [ "$idListFile" != "" ]; then 
    if [ -s "$idListFile" ]; then 
        idList="$idList `/bin/cat $idListFile`"
    else
        echo Warning! idListFile \'$idListFile\' does not exist or empty. >&2
    fi
fi

if [ ${#idList} -eq 0 ]; then
    echo No pfamid set. Exit. >&2
    exit 1
fi

if [ "$outpath" == "" ]; then 
    echo Outpath not set. Exit >&2
    exit 1
fi

mkdir -p $outpath

echo Command line: $commandline
echo Start at    `/bin/date`
general_res1=$(/bin/date +%s.%N)

if [ ! -d "$datapath"  ] ; then
    if [ "$dbnameseq" == "" ]; then 
        echo Neither datapath nor dbnameseq is set. Exit >&2
        exit 1
    elif [ ! -f "$dbnameseq.indexbin" -o ! -f "$dbnameseq.index" ]; then
        echo dbfile for \"$dbnameseq\" does not exist. Exist >&2
        exit 1
    fi
    
    tmpdir=$(mktemp -d /tmp/tmpdir.$progname.XXXXXXXXX) \
        || { echo Failed to create temp dir >&2; exit 1; }  
    datapath=$tmpdir
    tmpidlistfile=$tmpdir/tmpidlist.txt
    echo $idList | tr ' ' '\n' > $tmpidlistfile

    case $anamode in 
        0)
            if [ "$dbnameseq" == "" ]; then
                echo dbnameseq. Exit. >&2
                exit
            fi
            python $binpath/my_extractdb.py -dbname $dbnameseq -l $tmpidlistfile \
                -split -dataext .fa -outpath $tmpdir
            ;;
        *)  
            echo Unrecognized anamode = $anamode. Exit >&2
            exit
            ;;
    esac
fi

((cnt=0))
for pfamid in $idList ; do 
    if [ $cnt -ge $begin -a $cnt -lt $end ]; then 
# add a layer to show the running time
        res1=$(/bin/date +%s.%N)
        case $anamode in 
            0) SelectPairWithinPfam $pfamid;;
        esac
#        sleep 1s
        res2=$(/bin/date +%s.%N)
        printf "%s: Running time for %s:    %.3F\n" "$progname" \
            $pfamid  $(echo "$res2 - $res1"|/usr/bin/bc )
        echo
    fi
    ((cnt++))
done                                                                 
echo Finished at `/bin/date`
general_res2=$(/bin/date +%s.%N)
printf "%s: Running time for %d items is: %.3F\n" "$progname" $cnt  \
    $(echo "$general_res2 - $general_res1"|/usr/bin/bc )

# dump together
cd $outpath
cat PF*tableinfo.txt > allpfam.tableinfo.txt
cat PF*pairaln.fa > allpfam.pairaln.fa
rm -f $errFile
rm -rf $tmpdir

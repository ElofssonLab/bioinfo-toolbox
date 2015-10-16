#!/bin/bash
# Analyzing the topology variation of alpha helical transmembrane proteins
# This script is is a continued version of the script pfamTopoAna.sh since
# 2011-11-19 

# ChangeLog 2011-11-16#{{{
#    1. Cleaning should be done before multiple sequence alignment 
#    2. If there are sequences dropped in cleaning, if the input is a
#       multiple sequence alignment, unnecessarily gaps should be removed.
# ChangeLog 2011-11-19 
#    1. Bug solved in AnaMSATopo2, the following command is run when
#    randfasta.py is not run
#       /bin/cp -f $tmpnrfastafile $limitedFastaFile_cleaned
# ChangeLog 2012-02-23 
#    1. For mode2, alignment programs can be used as option
# ChangeLog 2012-03-08
#    when running python script, add python before that, so that the default
#    python module will be loaded.
# ChangeLog 2012-04-12
#    option topopredprog added. topology can also be predicted by
#    topcons_single
# ChangeLog 2012-11-08
#    The name of the topocons_single predicted topology output file fixed
# ChangeLog 2013-04-16
#    renamed_msaInFastaFormat is cleaned after running
# ChangeLog 2014-09-26
#    add the option, -forcewrite
#}}}

progname=`basename $0`
size_progname=${#progname}
wspace=`printf "%*s" $size_progname ""` 

rundir=`dirname $0`
anamode=2  #default anamode = 2
binpath=$rundir
dgpredpath=$rundir/../../../program/dgpred_standalone
newscampiscriptpath=$rundir/../../../program/newscampiscript
run_kalignP_path=$rundir/../../../program/run_kalignP scampi_dir=$rundir/../../../share/scampi
modhmm_bin=$rundir/../../../share/modhmm/bin
emboss_bin=$rundir/../../../usr/share/emboss/bin
cdhit_bin=$rundir/../../../bin
clustalo_bin=$rundir/../../../bin
kalign_bin=$rundir/../../../bin
fasttree_bin=$rundir/../../../bin
isPrintVerboseInfo= # in bash, if [ $istrue ]. the empty variable means false and any non-empty variable means true 
blastbin=$BLASTBIN
#blastdb=$BLASTDB/uniprotWG100Filtered
blastdb=$BLASTDB/uniprotWG100
#blastdb=$BLASTDB/swissprot


MAX_NUM_SEQ=1000 # maximum number of sequences to be processed in the multiple sequence alignment for anamode3.


usage="
Usage: $progname [-b INT] [-e INT] [-v] [-outpath]
       $wspace [-anamode 0   -dbnameseq STR  -dbnamemsa STR  ]
       $wspace [-anamode 1|2 -dbnameseq STR ]
       $wspace [-anamode 3   -dbnameseq STR  [ -dbnametopo STR] ]
       $wspace [-datapath DIR]
       $wspace -l LISTFILE  or ID [ID...]

Options:
    -l FILE           Set the file containing ids (rootname) of sequences
    -b INT            Begin index for id, default = 0
    -e INT            End index for id, default = INF
    -outpath  DIR     Location to output the result, (default: ./)
    -anamode  INT     Set the analysis mode, 0, 1, 2 or 3 (default: $anamode)
                      0. given sequences and multiple alignments
                      1. given sequences, alignment is obtained by kalign
                      2. given sequences, alignment is kalignP
                      3. given query sequence, homologous is obtained by blast
                         searching
    -alnprog STR      Set alignment program, (default: kalignp)
                      can be kalign, kalignp, clustalo
    -topopredprog STR Topology prediction method, (default: topocons_single)
                      can be scampi, topcons_single
    -noblast          Do not run psibalst. However, if \$id.blast does not exist.
                      Blast will force to run.
    -dbnameseq STR    Set the database name for sequences
    -dbnamemsa STR    Set the database name for MSAs in msf format
    -dbnametopo STR   Set the database name for topology
    -datapath  DIR    Set the datapath, if supplied, 
                      Sequence files will be searched at \$datapath/\$id.fa,
                      Topology file at \$datapath/\$id.topo
                      Msafile at \$datapath/\$id.msa.fa
    -startfrom INT    Start from step N.
    -forcewrite       Force overwrite the existing result
    -min-rlty  FLOAT  Miminum reliability score
    -mcmp      INT    method for compareMSATopo
    -maxseq    INT    Set maximum number of sequences to be kept in the
                      multiple alignment, (default: 1000)
    -extfa    STR     File extension for sequence file, (default: .fa)
    -exttopo  STR     File extension for topology file, (default: .topo)
    -extmsa   STR     File extension for multiple sequence alignment file, (default: .msa.fa)
    -v, --verbose     Print verbose information
    -h, --help        Print this help message and exit

Created 2010-08-17, updated 2014-09-26, Nanjiang Shu
nanjiang.shu@scilifelab.se
"
#export  GDFONTPATH=$DATADIR3/fonts/msttcorefonts/ 

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
exec_cmd(){ #{{{
    echo "$*"
    eval "$*"
}
#}}}
RunTopoAna(){ #{{{
    local id=$1
    # add a layer to show the running time
    res1=$(/bin/date +%s.%N)
    case $anamode in 
        0) AnaMSATopo0 $id;;
        1) AnaMSATopo1 $id;;
        2) AnaMSATopo2 $id;;
        3) AnaMSATopo3 $id;;
    esac
    #        sleep 1s
    res2=$(/bin/date +%s.%N)
    printf "%s: Running time for %s:    %.3F\n" "$progname" \
        $id  $(echo "$res2 - $res1"|/usr/bin/bc )
    echo
}
#}}}
AnaMSATopo0(){ # $id#{{{ # using multiple sequence alignment from pfam,  sequences are incomplete
    local id=$1
 
    local fastaFile=$datapath/$id.fa
    local msaFile=$datapath/$id.msf
    local msaInFastaFormat=$outpath/$id.msa.fa
    if [ ! -f "$msaFile" ] ; then
        echo msaFile \'$msaFile\' does not exist for $id. Ignore. >&2
        return 1
    fi
    $emboss_bin/seqret -sf msf -sequence $msaFile -osf fasta -outseq $msaInFastaFormat

    ### Step 1: MP topology prediction for each sequence 
    ### using scampi using the scampi without xslt
    $newscampiscriptpath/mySCAMPI_run.pl $fastaFile --scampipath $scampi_dir \
        --modhmmpath $modhmm_bin --outpath $outpath
    local topoFile=$outpath/$id.fa.topo
  
    ### Step 2: get MSATopoSeq
    local msatopoSeqFile=$outpath/$id.topomsa.fa
    $binpath/matchMSAtopo -msa $msaInFastaFormat  -topo $topoFile \
        -o $msatopoSeqFile
 
    # step 3: cleaning, remove single spanning and non TM proteins
    local msatopoSeqFile_cleaned=$outpath/$id.cleaned.topomsa.fa
    python $binpath/cleanSingleSpanTMPro.py $msatopoSeqFile -o $msatopoSeqFile_cleaned
    local topoFile_cleaned=$outpath/$id.cleaned.topo
    python $binpath/cleanSingleSpanTMPro.py $topoFile -o $topoFile_cleaned
 
    ### Step 3: plotting
    echo $binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile
    #$binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile

    echo not finished
    echo 
}
#}}}
AnaMSATopo1(){ # $id#{{{ # using kalign to get MSA
    #Note 2011-11-16
    # 1. Cleaning should be done before multiple sequence alignment 
    # 2. If there are sequences dropped in cleaning, if the input is a
    # multiple sequence alignment, unnecessarily gaps should be removed.
    local id=$1
    local fastaFile=$datapath/$id.fa
    if [ ! -s "$fastaFile" ]; then 
        echo seqfile \'$fastaFile\' does not exist or empty. Ignore. >&2
    else
        local numseq=`python $binpath/countseq.py $fastaFile -nf `
        if [ $numseq -le 1 ]; then 
            echo Too few sequences \($numseq\) for ID $id. \
                At least 2 sequences are needed. Ignore. >&2
            return 1
        fi
    fi

    local repeatedIDList=`IsAllSeqIDUnique $fastaFile`
    if [ "$repeatedIDList" != "" ] ; then
        echo -e Repeated sequence id found in file \'$fastaFile\'. Ignore. \
            They are\n"$repeatedIDList" >&2
        return 1
    fi
    ### Step 1: MP topology prediction for each sequence 
    ### using scampi  without xslt
    $newscampiscriptpath/mySCAMPI_run.pl $fastaFile --scampipath $scampi_dir \
        --modhmmpath $modhmm_bin --outpath $outpath
    local topoFile=$outpath/$id.fa.topo
    if [ ! -s "$topoFile" ] ; then
        echo Failed to predict topology for ID $id. Ignore. >&2
        return 1
    fi
    ### Step 2: cleaning, remove single spanning and non TM proteins
    local topoFile_cleaned=$outpath/$id.cleaned.topo
    python $binpath/cleanSingleSpanTMPro.py $topoFile -o $topoFile_cleaned
    if [ ! -s "$topoFile_cleaned" ]; then
        echo Failed to clean the topology for ID $id. Ignore. >&2
        return 1
    fi

    local tmpcleanedidlistfile=$outpath/$id.cleaned.idlist 
    local fastaFile_cleaned=$outpath/$id.cleaned.fa
    python $binpath/getfastaid.py $topoFile_cleaned -o $tmpcleanedidlistfile
    python $binpath/selectfastaseq.py -f $fastaFile -l $tmpcleanedidlistfile \
        -o $fastaFile_cleaned
    rm -f $tmpcleanedidlistfile

    local msaInFastaFormat=$outpath/$id.cleaned.msa.fa
    $BINPATH/kalign -i $fastaFile_cleaned -f fasta -o $msaInFastaFormat -q
  
    ### Step 2: get MSATopoSeq
    local msatopoSeqFile=$outpath/$id.cleaned.topomsa.fa
    $binpath/matchMSAtopo -msa $msaInFastaFormat -topo $topoFile \
        -o $msatopoSeqFile
 
    ### Step 5: compare topologies 
    local resultdifffile=$outpath/$id.diff.ana
    local sortedOrigTopoMSAFile=$outpath/$id.sorted.orig.topomsa.fa
    python $binpath/compareMSATopo.py $msaTopoSeqFile -o $resultdifffile \
        -wo $sortedOrigTopoMSAFile
}
#}}}
AnaMSATopo2(){ # $id#{{{  #using kalignP and a fasta file with multiple homologous sequences
    # Created 2011-08-31
    # Derived from AnaMSATopo1. run_kalignP is used instead of kalign.
    # Single-spanning TM proteins not not filtered since they may be caused by
    # deletion of a TM helix or by mis-prediction.
    # == updated 2011-11-16
    local id=$1
    local fastaFile=$datapath/${id}${ext_fa}

    if [ ! -s "$fastaFile" ]; then 
        echo seqfile \'$fastaFile\' does not exist or empty. Ignore. >&2
    else
        local numseq=`python $binpath/countseq.py $fastaFile -nf `
        if [ $numseq -le 1 ]; then 
            echo Too few sequences \($numseq\) for ID $id. At least 2 \
                sequences are needed. Ignore. >&2
            return 1
        fi
    fi

    local repeatedIDList=`IsAllSeqIDUnique $fastaFile`
    if [ "$repeatedIDList" != "" ] ; then
        echo -e Repeated sequence id found in file $fastaFile. Ignore. \
            They are\n"$repeatedIDList" >&2
        return 1
    fi

    ### Step 1: MP topology prediction for each sequence 
    ### using scampi using the scampi without xslt
    local topoFile=$datapath/${id}${ext_topo}
    if [ $startfrom -le 1 ]; then
        echo  Step 1: Get topologies for file $fastaFile...
        #if [ ! -s "$topoFile"  ]; then
        if [ ! -f "$topoFile"  ]; then
            echo "topoFile '$topoFile' does not exist, try to obtain it"
            if [  "$dbnametopo" != "" ];then
                echo "obtain topofile from the database $dbnametopo"
                local tmp_seqidfile=$outpath/.$id.tmp.seqidfile
                python $binpath/getfastaid.py -i $fastaFile -o $tmp_seqidfile
                python $binpath/my_extractdb.py -l $tmp_seqidfile \
                    -dbname $dbnametopo > $topoFile
                rm -f $tmp_seqidfile
            else
                echo "predict topology using $topologyPredProg"
                case $topologyPredProg in 
                    scampi|scampi_single)
                        $newscampiscriptpath/mySCAMPI_run.pl $fastaFile \
                            --scampipath $scampi_dir --modhmmpath $modhmm_bin \
                            --outpath $outpath
                        /bin/rm -f $outpath/$id.*.res
                        ;;
                    topcons_single)
                        local allinfofile=$outpath/$id.topcons-single.allinfo
                        $binpath/run_topcons_single.sh $fastaFile -outpath $outpath
                        python $binpath/topcons_single_to_fasta.py $allinfofile -outpath $outpath
                        mv $outpath/$id.topcons-single_topcons_single.topo $topoFile
                esac
            fi
        fi
    fi

    local msaFile=${datapath}/${id}${ext_msa}
    if [ ! -s "$msaFile" ];then
        ### Step 2: cleaning, remove non TM proteins
        local topoFile_cleaned=$outpath/$id.cleaned.topo
        if [ $startfrom -le 2 ]; then
            if [ ! -s "$topoFile" ] ; then
                echo Failed to predict topology for ID $id. Ignore. >&2
                return 1
            fi
            echo Step 2: Remove non TM proteins...
            exec_cmd "python $binpath/cleanSingleSpanTMPro.py $topoFile -mintm 1 -o $topoFile_cleaned"

            local tmp_topoFile=$topoFile.tmp
            /bin/cp -f $topoFile_cleaned $tmp_topoFile
            echo Step 3: remove TMpros predicted with reliability score smaller $MIN_RLTY
            exec_cmd "python $binpath/cleanTMPro_by_RLTY.py $tmp_topoFile -min-rlty $MIN_RLTY  -o $topoFile_cleaned"
            rm -f $tmp_topoFile
        fi
        # step 3. Get cleaned amino acid sequences
        local tmpcleanedidlistfile=$outpath/$id.cleaned.idlist 
        local fastaFile_cleaned=$outpath/$id.cleaned.fa
        if [ $startfrom -le 3 ]; then
            if [ ! -s "$topoFile_cleaned" ]; then
                echo Failed to clean the topology for ID $id. Ignore. >&2
                return 1
            fi
            python $binpath/getfastaid.py $topoFile_cleaned -o $tmpcleanedidlistfile
            python $binpath/selectfastaseq.py -f $fastaFile -l $tmpcleanedidlistfile \
                -o $fastaFile_cleaned
            rm -f $tmpcleanedidlistfile
        fi
        # Step 4 limit the number of sequences to MAX_NUM_SEQ
        local limitedFastaFile_cleaned=$outpath/$id.homology.cleaned.le${MAX_NUM_SEQ}.fa 
        if [ $startfrom -le 4 ]; then
            if [ ! -s $fastaFile_cleaned ] ; then 
                echo Failed to get cleaned fasta seq file for ID $id. Ignore. >&2
                return 1
            fi
            local tmpnrfastafile=$outpath/$id.homology.cdhit.nr.fa.tmp
            $cdhit_bin/cd-hit -i $fastaFile_cleaned -o $tmpnrfastafile -c 1.0
            if [ ! -s "$tmpnrfastafile" ]; then 
                echo cd-hit failed for ID $id. just copy the fastaFile_cleaned >&2
                /bin/cp -f $fastaFile_cleaned $tmpnrfastafile
            fi
            local numseq=`python $binpath/countseq.py -nf $tmpnrfastafile`
            if [ $numseq -le 1 ]; then 
                echo Too few \($numseq\) sequences in the cleaned topology file \
                    for ID $id. Ignore. >&2
                return 1
            else
                if [ $numseq -gt $MAX_NUM_SEQ ] ; then 
                    # decrease the sequence identity threshold until 
                    # threshold >=0.75 or numseq <= MAX_NUM_SEQ
                    local tmpfile=$(mktemp /tmp/tmp.$progname.XXXXXXXXX) \
                        || { echo Failed to create temp file >&2; exit 1; }   
                    min_pid_threshold=75
                    pid_threshold=95
                    while [ 1 ] ; do
                        if [ $pid_threshold -lt $min_pid_threshold ]; then 
                            break
                        fi
                        if [ $numseq -le $MAX_NUM_SEQ ]; then 
                            break
                        fi
                        seqidt_threshold=`$binpath/e $pid_threshold / 100.0` 
                        $cdhit_bin/cd-hit -i $tmpnrfastafile -o $tmpfile \
                            -c $seqidt_threshold
                        /bin/mv -f $tmpfile $tmpnrfastafile
                        ((pid_threshold-=5))
                        numseq=`python $binpath/countseq.py -nf $tmpnrfastafile`
                    done 

                    # if still too many sequences, randomly select $MAX_NUM_SEQ of
                    # sequences
                    if [ $numseq -gt $MAX_NUM_SEQ ]; then 
                        python $binpath/randfasta.py -i $tmpnrfastafile -n $MAX_NUM_SEQ \
                            -o $limitedFastaFile_cleaned
                    else
                        /bin/cp -f $tmpnrfastafile $limitedFastaFile_cleaned
                    fi
                    /bin/rm -f $tmpfile
                else 
                    /bin/cp -f $tmpnrfastafile $limitedFastaFile_cleaned
                fi
            fi
            #remove tmpfiles produced by cd-hit
            /bin/rm -f ${tmpnrfastafile}*
        fi

        local msaInFastaFormat=${limitedFastaFile_cleaned%.*}.kalignP.fasta
        if [ $startfrom -le 5 ]; then
            if [ ! -s "$limitedFastaFile_cleaned" ]; then
                echo Failed to get limited sequences for ID $id. Ignore. >&2
                return 1
            fi
            echo  Step 5: Running multiple sequence alignment for \
                $limitedFastaFile_cleaned...
            case $msaProg in 
                kalignp)
                    $run_kalignP_path/run_kalignP.sh  $limitedFastaFile_cleaned  \
                        -f fasta -outpath $outpath -q 2> $errFile
                    if [ -s $errFile -a ! -s $msaInFastaFormat ]; then
                        echo Run KalignP with PSGP failed for $id. Try with no PSGP >&2
                        $run_kalignP_path/run_kalignP.sh $limitedFastaFile_cleaned \
                            -no-psgp \
                            -f fasta -outpath $outpath -q 2> $errFile
                    fi
                    ;;
                clustalo)
                    $clustalo_bin/clustalo -i $limitedFastaFile_cleaned \
                        -o $msaInFastaFormat
                    ;;
                *)
                    $kalign_bin/kalign -f fasta $limitedFastaFile_cleaned \
                        -o $msaInFastaFormat
            esac
        fi
    else
        topoFile_cleaned=$topoFile
        local seqidListFile_TM=$outpath/${id}.TM.seqidlist
        $binpath/getfastaid.py $topoFile_cleaned -o $seqidListFile_TM
        local fastaFile_TM=$outpath/${id}.TM.fa
        local msaFile_TM=$outpath/${id}.TM.mfa
        exec_cmd "$binpath/selectfastaseq.py -f $fastaFile -l $seqidListFile_TM -o $fastaFile_TM"
        exec_cmd "$binpath/selectfastaseq.py -f $msaFile -l $seqidListFile_TM -o $msaFile_TM.tmp"
        exec_cmd "$binpath/removeUnnecessaryGap $msaFile_TM.tmp -o $msaFile_TM"

        fastaFile_cleaned=$fastaFile_TM
        msaInFastaFormat=$msaFile_TM
    fi

    ### Step 6: get MSATopoSeq
    local msatopoSeqFile=$outpath/$id.cleaned.topomsa.fa
    if [ $startfrom -le 6 ]; then
        if [ ! -s $msaInFastaFormat ] ; then 
            echo Failed to generate multiple sequence alignment for $id. Ignore. >&2
            return 1
        fi
        echo Step 6: Matching Sequence MSA to topology MSA...
        exec_cmd "$binpath/matchMSAtopo -msa $msaInFastaFormat -topo $topoFile_cleaned  -o $msatopoSeqFile"
    fi

    ### Step 7 create dg files
    if [ $startfrom -le 7 ]; then
        echo Step 7: create dg files
        exec_cmd "$binpath/getDGvalueTMOfTopo.sh -topo $msatopoSeqFile -aa $fastaFile_cleaned -outpath $outpath "
        /bin/rm -f $outpath/$id.*.dgscorelist
        rm -f $msatopoSeqFile
    fi


    ### Step 8: compare topologies 
    local topowithDGscoreFile=$outpath/$id.cleaned.topomsa.topowithdgscore
    local dgscoreFile=$outpath/$id.cleaned.topomsa.dgscore
    local sortedOrigTopoMSAFile=$outpath/$id.sorted.orig.topomsa.fa
    local clusteredOrigTopoMSAFile=$outpath/$id.clustered.orig.topomsa.fa
    local clusteredOrigTopoMSAAnnoFile=$outpath/$id.clustered.orig.topomsa.anno
    local resultfile=$outpath/$id.diff.ana
    local groupedSortedOrigTopoMSAFile=$outpath/$id.grouped.sorted.orig.topomsa.fa
    local groupedResultFile=$outpath/$id.grouped.diff.ana
    local invertedTopologyFile=$outpath/$id.inverted.info.txt
    maxdgdiff=1.0
    if [ $startfrom -le 8 ]; then
        echo Step 8: compare MSA topologies
        if [ ! -s $topowithDGscoreFile ] ; then 
            echo Failed to generate topowithDGscoreFile for $id. Ignore. >&2
            return 1
        fi
        exec_cmd "python $binpath/compareMSATopo.py $topowithDGscoreFile -mcmp $method_comparison \
            -maxdgdiff $maxdgdiff \
            -wo $sortedOrigTopoMSAFile -o $resultfile \
            -og $groupedResultFile -wog $groupedSortedOrigTopoMSAFile \
            -woc $clusteredOrigTopoMSAFile \
            -woinv $invertedTopologyFile
            "
        grep "^>" $clusteredOrigTopoMSAFile > $clusteredOrigTopoMSAAnnoFile
        rm -f $clusteredOrigTopoMSAFile

        if [ -s "$topowithDGscoreFile" ]; then 
            awk '{if(NR%3== 0 || NR%3 == 1) print}' $topowithDGscoreFile > $dgscoreFile
        fi
        rm -f $topowithDGscoreFile
    fi

    #return 0
    # calculate KR bias 
    local invertedTopologyIDListFile=$invertedTopologyFile.INVseq.idlist
    local invertedTopologyFastaSeqFile=$invertedTopologyFile.INVseq.fa
    if [ -s "$invertedTopologyIDListFile" ];then
        exec_cmd "$binpath/selectfastaseq.py -l $invertedTopologyIDListFile -f $fastaFile -o $invertedTopologyFastaSeqFile"
        exec_cmd "$binpath/test_calculate_KR_bias.py -maxdist 12 -flankwin 0 -seqfile $invertedTopologyFastaSeqFile -topofile $topoFile -outpath $outpath"
        local krbiasfile=${invertedTopologyFastaSeqFile%.*}.krbias.txt
        if [ -s "$krbiasfile" ];then
            exec_cmd "$binpath/plotKRbiasHistogram.sh $krbiasfile"
        else
            echo "krbiasfile $krbiasfile does not exist"
        fi
    fi

    ### Step 9: Create pictures
    # create thumbnail from non text orig
    if [ $startfrom -le 9 ]; then
        echo Step 9: Draw figures
        if [ ! -s $sortedOrigTopoMSAFile -a ! -s $groupedSortedOrigTopoMSAFile ]; then
            echo Failed to generate sorted TopoMSA file for $id. >&2
            return 1
        fi
        exec_cmd "python $binpath/drawMSATopo.py -text n -outpath $outpath\
            $sortedOrigTopoMSAFile \
            $groupedSortedOrigTopoMSAFile -aapath $datapath -showTMidx -colorTMbox -pfm no"

        /usr/bin/convert -thumbnail 200 $outpath/$id.sorted.orig.topomsa.png \
            $outpath/thumb.$id.sorted.orig.topomsa.png
        /usr/bin/convert -thumbnail 200 $outpath/$id.grouped.sorted.orig.topomsa.png \
            $outpath/thumb.$id.grouped.sorted.orig.topomsa.png

        # create full size image with text (text is amino acid sequence)
        exec_cmd "python $binpath/drawMSATopo.py -text y -outpath $outpath\
            $sortedOrigTopoMSAFile \
            $groupedSortedOrigTopoMSAFile\
            -aapath $datapath -showTMidx -colorTMbox -pfm no"

        #rm -f $topoFile
        #rm -f $topoFile_cleaned
        #rm -f $fastaFile_cleaned
        #rm -f $limitedFastaFile_cleaned
    fi

    # create tree
    local renamed_msaInFastaFormat=$outpath/$id.renamedid.msa.fasta
    $binpath/renameSeqIDInFasta.py $msaInFastaFormat -o $renamed_msaInFastaFormat
    treeFile=$outpath/$id.kalignp.fasttree
    if [  ! -s $treeFile  -o "$isOverwrite" == "1" ]; then 
        exec_cmd "$fasttree_bin/FastTree $renamed_msaInFastaFormat > $treeFile"
    fi
    $binpath/sortedTopoMSA2colordef.sh $sortedOrigTopoMSAFile > $outpath/$id.cmpclass.colordef.txt
    $binpath/clusteredTopoMSA2colordef.sh $clusteredOrigTopoMSAAnnoFile > $outpath/$id.cluster.colordef.txt
    $binpath/sortedTopoMSA2numTMbardef.sh $sortedOrigTopoMSAFile > $outpath/$id.numTMdef.txt
    python $binpath/sortedTopoMSA2inside-outside-colordef.py $sortedOrigTopoMSAFile > $outpath/$id.ntermstate.colordef.txt
    python $binpath/sortedTopoMSA2numTM_and_io.py $sortedOrigTopoMSAFile > $outpath/$id.numTM_and_io.txt
    python $binpath/itol_pfamtree.py -datapath $outpath -outpath $outpath $id
    rm -f $renamed_msaInFastaFormat

    # draw reordered 
    local reorderedmsafile=$outpath/${id}.reordered.topomsa.fa
    local bsname_treeFile=`basename $treeFile`
    local rtname_treeFile=${bsname_treeFile%.*}
    local orderlistfile=$outpath/${rtname_treeFile}.listorder.txt
    exec_cmd "$binpath/itol_get_tree_listorder.py $treeFile -outpath $outpath"
    if [ ! -s "$sortedOrigTopoMSAFile" ]; then
        echo "msafile $sortedOrigTopoMSAFile does not exist or empty. Ignore $id" >&2
        return 1
    fi
    if [ ! -s "$orderlistfile" ]; then
        echo "orderlistfile $orderlistfile does not exist or empty. Ignore $id" >&2
        return 1
    fi
    exec_cmd "python $binpath/reordermsa.py -msafile $sortedOrigTopoMSAFile -orderlist $orderlistfile -o $reorderedmsafile"
    if [ -f "$reorderedmsafile" ]; then
        exec_cmd "python $binpath/drawMSATopo.py -sep n -text n -aapath $datapath  -showTMidx -colorTMbox -pfm no -outpath $outpath $reorderedmsafile"
        exec_cmd "/usr/bin/convert -thumbnail 200 $outpath/${id}.reordered.topomsa.png $outpath/thumb.${id}.reordered.topomsa.png"
        #rm -f $reorderedmsafile
    fi

    echo
}
#}}}
AnaMSATopo3(){ #{{{ #$id #using kalignP and sequence and topology of query is given
    # Created 2011-09-15, updated 2011-09-15, Nanjiang Shu 
    # using query sequence and blast search to get homologous
    # Single-spanning TM proteins are not filted since they may be caused by
    # deletion of a TM helix, or by mis-prediction.
    # The topology of the query sequencce will be predicted my Scampi if not
    # given 
    # == updated 2011-11-16
    local id=$1
    echo -e Analyzing ID $id\n

    #Step 1: get query sequence and topology files
    echo -e Step 1: Get query sequence and query topology...\n
    local querySeqFile=$datapath/$id.fa
    if [ ! -s "$querySeqFile" ] ; then 
        echo query sequence file \"$querySeqFile\" does not exist. Ignore. >&2
        return 1
    fi
    local queryTopoFile=$datapath/$id.topo
    if [ ! -s "$queryTopoFile" ] ; then 
        $newscampiscriptpath/mySCAMPI_run.pl $querySeqFile \
            --scampipath $scampi_dir --modhmmpath $modhmm_bin --outpath $outpath
        mv -f $outpath/$id.fa.topo $outpath/$id.topo
        rm -f $outpath/$id.*.res
    else
        /bin/cp -f $queryTopoFile $outpath
    fi
    local queryTopoFile=$outpath/$id.topo

    ### Step 2: obtain the fasta sequences for those with E-values <= the
    ### threshold
    ### obtain only hits from the last round, for multiple hits from the same
    ### sequence, take the one with lowest evalue
    ### sequences are obtained by fastacmd
    local blastfile=$outpath/$id.blast
    echo -e Step 2: Get homologous sequences using psiblast...\n
    if [ $isRunBlast -eq 1 -o ! -s $blastfile ]; then 
        $blastbin/blastpgp -i $querySeqFile -h 1e-5 -d $blastdb -j 3 -a 1 \
            -o $blastfile -m 9 -b 10000
    fi 

    local homologyFastaFile=$outpath/$id.origblast.fa
    local evalue=1e-5
    python $binpath/blastm9tofasta.py $blastfile -evalue $evalue -blastdb $blastdb \
        -o $homologyFastaFile

    ### Step 3: Membrane protein topology prediction for each sequence in homologyFastaFile
    ### using scampi using the scampi without xslt
    echo -e Step 3: Predict membrane protein topology by Scampi for \
        file $homologyFastaFile...\n
    $newscampiscriptpath/mySCAMPI_run.pl $homologyFastaFile \
        --scampipath $scampi_dir --modhmmpath $modhmm_bin --outpath $outpath
    rm -f $outpath/$id.*.res
    local homologyTopoFile=$homologyFastaFile.topo

    if [ ! -s "$homologyTopoFile" ]; then 
        echo Failed to get topology file for homologous sequences for ID $id.\
            Ignore. >&2
        return 1
    fi

    local topoFile=$outpath/$id.homology.topo
    cat $queryTopoFile $homologyTopoFile > $topoFile

    local fastaFile=$outpath/$id.homology.fa
    cat $querySeqFile $homologyFastaFile > $fastaFile

    ### Step 4: cleaning non TM proteins
    echo -e Step 4: Clean non TM proteins...\n
    local topoFile_cleaned=$outpath/$id.homology.cleaned.topo
    python $binpath/cleanSingleSpanTMPro.py $topoFile -mintm 1 -o $topoFile_cleaned
    if [ ! -s "$topoFile_cleaned" ]; then 
        echo Failed to get filtered topology MSA file for ID $id. Ignore. >&2
        return 1
    fi

    # step 4.1: Get cleaned amino acid sequences
    local tmpcleanedidlistfile=$outpath/$id.homology.cleaned.idlist 
    local fastaFile_cleaned=$outpath/$id.homology.cleaned.fa
    python $binpath/getfastaid.py $topoFile_cleaned -o $tmpcleanedidlistfile
    python $binpath/selectfastaseq.py -f $fastaFile -l $tmpcleanedidlistfile -o $fastaFile_cleaned
    rm -f $tmpcleanedidlistfile
    if [ ! -s $fastaFile_cleaned ] ; then 
        echo Failed to get cleaned fasta seq file for ID $id. Ignore. >&2
        return 1
    fi

    # Step 5: limit the number of sequences to MAX_NUM_SEQ
    echo -e Step 5: Limit the number of sequences to $MAX_NUM_SEQ...\n
    #5.1 cd-hit to remove redundant sequences
    local tmpnrfastafile=$outpath/$id.homology.cdhit.nr.fa.tmp
    $cdhit_bin/cd-hit -i $fastaFile_cleaned -o $tmpnrfastafile -c 1.0
    if [ ! -s $tmpnrfastafile ]; then 
        echo cd-hit failed. just copy the fastaFile_cleaned >&2
        /bin/cp -f $fastaFile_cleaned $tmpnrfastafile
    fi
    local numseq=`python $binpath/countseq.py -nf $tmpnrfastafile`
    if [ $numseq -le 1 ]; then 
        echo Too few \($numseq\) sequences in the cleaned topology file \
            for ID $id. Ignore. >&2
        return 1
    else
        local limitedFastaFile_cleaned=$outpath/$id.homology.cleaned.le${MAX_NUM_SEQ}.fa 
        if [ $numseq -gt $MAX_NUM_SEQ ] ; then 
            # decrease the sequence identity threshold until 
            # threshold >=0.75 or numseq <= MAX_NUM_SEQ
            local tmpfile=$(mktemp /tmp/tmp.$progname.XXXXXXXXX) \
                || { echo Failed to create temp file >&2; exit 1; }   
            min_pid_threshold=75
            pid_threshold=95
            while [ 1 ] ; do
                if [ $pid_threshold -lt $min_pid_threshold ]; then 
                    break
                fi
                if [ $numseq -le $MAX_NUM_SEQ ]; then 
                    break
                fi
                seqidt_threshold=`$binpath/e $pid_threshold / 100.0` 
                $cdhit_bin/cd-hit -i $tmpnrfastafile -o $tmpfile \
                    -c $seqidt_threshold
                /bin/mv -f $tmpfile $tmpnrfastafile
                ((pid_threshold-=5))
                numseq=`python $binpath/countseq.py -nf $tmpnrfastafile`
            done 

            # if still too many sequences, randomly select $MAX_NUM_SEQ of
            # sequences
            if [ $numseq -gt $MAX_NUM_SEQ ]; then 
                python $binpath/randfasta.py -i $tmpnrfastafile -n $MAX_NUM_SEQ \
                    -o $limitedFastaFile_cleaned
            else
                /bin/cp -f $tmpnrfastafile $limitedFastaFile_cleaned
            fi
            /bin/rm -f $tmpfile
        else 
            /bin/cp -f $tmpnrfastafile $limitedFastaFile_cleaned
        fi

        if [ ! -s "$limitedFastaFile_cleaned" ]; then
            echo Failed to get limited sequences for ID $id. Ignore. >&2
            return 1
        fi
    fi

    ### Step 6: Get multiple sequence alignment
    echo -e Step 6: Do multiple sequence alignment for \
        $limitedFastaFile_cleaned...\n
    local msaInFastaFormat=${limitedFastaFile_cleaned%.*}.kalignP.fasta
    $run_kalignP_path/run_kalignP.sh  $limitedFastaFile_cleaned -f fasta \
        -outpath $outpath -q
    
    if [ ! -s $msaInFastaFormat ]; then
        echo Run KalignP with PSGP failed for $fastaFile. Try with no PSGP >&2
        $run_kalignP_path/run_kalignP.sh  $limitedFastaFile_cleaned -no-psgp \
            -f fasta -outpath $outpath -q 2> $errFile
    fi

    if [ ! -s $msaInFastaFormat ] ; then 
        echo Failed to generate multiple sequence alignment \
            for $limitedFastaFile_cleaned. Ignore. >&2
        return 1
    fi
  
    ### Step 7: get MSATopoSeq
    local msatopoSeqFile=$outpath/$id.homology.cleaned.topomsa.fa
    echo -e Step 7: Get topology MSA by matching with sequence MSA...\n
    $binpath/matchMSAtopo -msa $msaInFastaFormat -topo $topoFile_cleaned \
        -o $msatopoSeqFile
    if [ ! -s "$msatopoSeqFile" ]; then 
        echo Failed to match topology MSA for ID $id. Ignore. >&2
        return 1
    fi
 
    ### Step 7 create dg files
    echo -e Step 7: Add DG values of TM helices to topology MSA...\n
    $binpath/getDGvalueTMOfTopo.sh -topo $msatopoSeqFile \
        -aa $limitedFastaFile_cleaned -outpath $outpath  
    /bin/rm -f $outpath/$id.*.dgscorelist
    local topowithDGscoreFile=$outpath/$id.homology.cleaned.topomsa.topowithdgscore

    if [ ! -s "$topowithDGscoreFile" ]; then 
        echo Failed to add DG score for ID $id. Ignore.
        return 1
    fi

    ### Step 8: compare topologies 
    echo -e Step 8: Compare topologies in the topology MSA file \
        $topoWithDGscoreFile...\n
#     local sortedOrigTopoMSAFile=$outpath/$id.homology.sorted.orig.topomsa.fa
#     local resultfile=$outpath/$id.homology.diff.ana
    local topowithDGscoreFile=$outpath/$id.homology.cleaned.topomsa.topowithdgscore
    local sortedOrigTopoMSAFile=$outpath/$id.homology.sorted.orig.topomsa.fa
    local clusteredOrigTopoMSAFile=$outpath/$id.homology.clustered.orig.topomsa.fa
    local clusteredOrigTopoMSAAnnoFile=$outpath/$id.homology.clustered.orig.topomsa.anno
    local resultfile=$outpath/$id.homology.diff.ana
    local groupedSortedOrigTopoMSAFile=$outpath/$id.homology.grouped.sorted.orig.topomsa.fa
    local groupedResultFile=$outpath/$id.homology.grouped.diff.ana
    maxdgdiff=1.0
    #python $binpath/compareMSATopo.py $topowithDGscoreFile  -mcmp $method_comparison -wo $sortedOrigTopoMSAFile  -o $resultfile 
    python $binpath/compareMSATopo.py $topowithDGscoreFile -mcmp $method_comparison \
        -maxdgdiff $maxdgdiff \
        -wo $sortedOrigTopoMSAFile -o $resultfile \
        -og $groupedResultFile -wog $groupedSortedOrigTopoMSAFile \
        -woc $clusteredOrigTopoMSAFile
    grep "^>" $clusteredOrigTopoMSAFile > $clusteredOrigTopoMSAAnnoFile
    rm -f $clusteredOrigTopoMSAFile

    ### Step 9: Draw figures
    # create thumbnail from non text orig
    echo -e Step 9: Draw figures for the topology MSA...\n
    python $binpath/drawMSATopo.py -text n $sortedOrigTopoMSAFile -outpath $outpath
    /usr/bin/convert -thumbnail 200 \
        $outpath/$id.homology.sorted.orig.topomsa.png \
        $outpath/thumb.$id.homology.sorted.orig.topomsa.png

    # create full size image with text (text is amino acid sequence)
    python $binpath/drawMSATopo.py -text y $sortedOrigTopoMSAFile -outpath $outpath
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
dbnameseq=
dbnamemsa=
dbnametopo=
datapath=
outpath=
begin=0
end=9999999
msaProg=kalignp
topologyPredProg=topcons_single


method_comparison=1
isDrawText=yes
idList=
isRunBlast=1
startfrom=0

MIN_RLTY=0

ext_fa=.fa
ext_topo=.topo
ext_msa=.msa.fa


isNonOptionArg=0
isOverwrite=0
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
            -dbnamemsa|--dbnamemsa) dbnamemsa=$2;shift;;
            -dbnametopo|--dbnametopo) dbnametopo=$2;shift;;
            -maxseq|--maxseq) MAX_NUM_SEQ=$2;shift;;
            -extfa|--extfa) ext_fa=$2;shift;;
            -exttopo|--exttopo) ext_topo=$2;shift;;
            -forcewrite|--forcewrite) isOverwrite=1;;
            -extmsa|--extmsa) ext_msa=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -noblast|--noblast) isRunBlast=0;;
            -alnprog|--alnprog) msaProg=$2;shift;;
            -topopredprog|--topopredprog) topologyPredProg=$2;shift;;
            -anamode|--anamode) anamode=$2;shift;;
            -min-rlty|--min-rlty) MIN_RLTY=$2;shift;;
            -mcmp|--mcmp) method_comparison=$2;shift;;
            -startfrom|--startfrom) startfrom=$2;shift;;
            -v|-verbose|--verbose) isPrintVerboseInfo=1;shift;;
            -*) echo Error! Wrong argument: $1; exit;;
        esac
    else
        idList="$idList $1"
    fi
    shift
done

if [ "$anamode" == "0" ];then
    IsProgExist $emboss_bin/seqret
fi
IsProgExist $newscampiscriptpath/mySCAMPI_run.pl
IsProgExist $BINPATH/kalign
IsProgExist /usr/bin/bc
IsProgExist /bin/date
IsProgExist $binpath/my_extractdb.py
IsProgExist $binpath/countseq.py
IsPathExist $run_kalignP_path
IsProgExist $BINPATH/e

if [ "$idListFile" != "" ]; then 
    if [ -s "$idListFile" ]; then 
        idList="$idList `/bin/cat $idListFile`"
    else
        echo Warning! idListFile \'$idListFile\' does not exist or empty. >&2
    fi
fi

if [ ${#idList} -eq 0 ]; then
    echo No ID set. Exit. >&2
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
    fi
    tmpdir=$(mktemp -d /tmp/tmpdir.$progname.XXXXXXXXX) \
        || { echo Failed to create temp dir >&2; exit 1; }  
    datapath=$tmpdir
    tmpidlistfile=$tmpdir/tmpidlist.txt
    echo $idList | tr ' ' '\n' > $tmpidlistfile

    case $anamode in 
        0)
            if [ "$dbnameseq" == "" -o "$dbnamemsa" == "" ]; then
                echo dbnameseq or dbnamemsa not set. Exit. >&2
                exit
            fi
            python $binpath/my_extractdb.py -dbname $dbnameseq -l $tmpidlistfile \
                -split -dataext .fa -outpath $tmpdir
            $binpath/my_extractdb.py -dbname $dbnamemsa -l $tmpidlistfile \
                -split -dataext .msf -outpath $tmpdir
            ;;
        1|2)
            if [ "$dbnameseq" == "" ]; then
                echo dbnameseq not set. Exit. >&2
                exit 1
            fi
            python $binpath/my_extractdb.py -dbname $dbnameseq -l $tmpidlistfile \
                -split -dataext .fa -outpath $tmpdir
            ;;
        3)
            if [ "$dbnameseq" == ""  ]; then
                echo dbnameseq not set. Exit. >&2
                exit 1
            fi
            python $binpath/my_extractdb.py -dbname $dbnameseq -l $tmpidlistfile \
                -split -dataext .fa -outpath $tmpdir
            if [ "$dbnametopo" != "" ]; then 
                python $binpath/my_extractdb.py -dbname $dbnametopo -l $tmpidlistfile \
                    -split -dataext .topo -outpath $tmpdir
            fi
            ;;
        *)  
            echo Unrecognized anamode = $anamode. Exit >&2
            exit
            ;;
    esac
fi

((cnt=0))
for id in $idList ; do 
    if [ $cnt -ge $begin -a $cnt -lt $end ]; then
        RunTopoAna $id
    fi
    ((cnt++))
done

echo Finished at `/bin/date`
general_res2=$(/bin/date +%s.%N)
printf "%s: Running time for %d items is: %.3F\n" "$progname" $cnt  \
    $(echo "$general_res2 - $general_res1"|/usr/bin/bc )

rm -f $errFile
rm -rf $tmpdir

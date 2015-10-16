#!/bin/bash
# analyzing membrane topology variation among homologous proteins
# 2010-08-17 

anamode=3  #default anamode = 1
binpath=/misc/illergarddata/wk/MPTopo/bin
dgpredpath=/misc/illergarddata/program/dgpred/dgpred_standalone
isPrintVerboseInfo= # in bash, if [ $istrue ]. the empty variable means false and any non-empty variable means true 

usage="
Usage: runMSATopoAna.sh -l idlistFile --aapath aaPath --blastpath blastPath
Options:
    -l          file : set the file containing ids (rootname) of sequences
    --aapath    path : set the path for amino acid sequence files, filename=\$aaPath/\$id.aa
    --blastpath path : set the path for blast files (output with m=0 option),filename=\$blastPath/\$id.blast
    --outpath   path : location to output the result
    --notuseorig     : do not use the orig sequence, but sequence in the local blast alignemnt
    --anamode 0|1|2|3: set the analysis mode, default=$anamode
                     : 0. run unlabelbed prediction and do topology comparisons
                     : 1. using scampi_msa to predict the query and run labeled prediction for homologues
    -v|--verbose     : print verbose information
    -h|--help        : print this help message and exit

Created 2010-08-17, updated 2010-09-13, Nanjiang
nanjiang.shu@gmail.com
"
export  GDFONTPATH=/data3/fonts/msttcorefonts/ 
function PrintHelp()
{
    echo "$usage"
}

function AnaMSATopo0() # $id#{{{
{
    local id=$1
    ### step 1: obtain the fasta sequences for those with E-values above the threshold
    blastFile=$blastPath/$id.blast
    ### extract the last round blast result
    lastRoundBlastFile=$outpath/$id.lastround.blast
    echo "extract_blast_round_n $blastFile -o  $lastRoundBlastFile "
    extract_blast_round_n $blastFile -o  $lastRoundBlastFile 
 
    tmpFastaFile=$outpath/$id.blast.tmp.fa
    blast2fasta.pl -e 1e-3 $lastRoundBlastFile > $tmpFastaFile
    
# remove duplicated IDs, this is caused by multiple HSPs from one sequence in
# the blast result
    runname=$outpath/$id.blast.tmp.uniq
    uniqueseq $runname -db $tmpFastaFile 
    uniqFastaFile=$outpath/$id.blast.tmp.uniq.seq

    tmpFastaFile=$uniqFastaFile
    blastedIDListFile=$outpath/$id.blast.tmp.idlist
    getfastaid.py -i $tmpFastaFile -o $blastedIDListFile
    if [ "$isUseOrigSeq" == "true" ]; then
        tmpFastaOrigFile=$outpath/$id.blast.tmp.fastacmd.fa
        fastacmd -i $blastedIDListFile -d uniprotWG100 1> $tmpFastaOrigFile 2> $errFile
        if [ -s $errFile ]; then 
            echo "fastacmd error, when running"
            echo "fastacmd -i $blastedIDListFile -d uniprotWG100 1> $tmpFastaOrigFile 2> $errFile"
            cat $errFile
        fi
        tmpAnnoFile=$outpath/$id.blast.tmp.anno
        $binpath/getfastaAnnotation.py -i $tmpFastaFile -o $tmpAnnoFile
        $binpath/replaceFastaAnnotation.py -fasta $tmpFastaOrigFile -anno $tmpAnnoFile -o $tmpFastaFile
    fi

    #in case there are too many sequences, reduce the number of sequences by culling at certain sequence identity level
    Nseq=`wc -l $blastedIDListFile | awk '{print $1}'`
    if [ $Nseq -ge 2000 ]; then
        tmpCDHITFile=$outpath/$id.blast.cdhit.tmp.fa
        cd-hit -i  $tmpFastaFile  -o $tmpCDHITFile -c 0.7 -n 4   
        tmpFastaFile=$tmpCDHITFile
        Nseq=`getfastaid.py -i $tmpFastaFile | wc -l`
        if [ $Nseq -ge 2000 ]; then
#            output only the first 2000
             top2000FastaFile=$outpath/$id.tmp.top2000.fa
             catfasta.py -b 0 -e 2000 $tmpFastaFile -o $top2000FastaFile 
             tmpFastaFile=$top2000FastaFile
        fi
    fi 

    echo "$id: fasta file from blast output generated!"
 
    fastaFile=$outpath/$id.blast.fa
    ### add the target sequence to the fasta file
    cat $aaPath/$id.aa $tmpFastaFile > $fastaFile
 
    ### Step 2: Kalign of these sequences
    msaFile=$outpath/$id.msa
    kalign -i $fastaFile -o $msaFile
 
    ### Step 3: MP topology prediction for each sequence 
    ### using scampi or dG_toppred  , using the scampi without xslt
    $CASIODATA3/bin/mySCAMPI_run.pl $fastaFile --outpath $outpath
    resFile=$outpath/$id.blast.fa.res
    topologyFile=$outpath/$id.blast.topo
#     echo "/scampi/compacttopo2top.py $resFile > $topologyFile"
#     /scampi/compacttopo2top.py $resFile > $topologyFile
  
    ### Step 4: get MSATopoSeq
    msatopoSeqFile=$outpath/$id.topomsa
    $binpath/matchMSAtopo -msa $msaFile -fa $fastaFile -topo $topologyFile > $msatopoSeqFile
    $binpath/matchMSAtopo -msa $msaFile -fa $fastaFile -topo $topologyFile --format 1 > $msatopoSeqFile.ver
 
    ### Step 5: pairwise comparison of topologies
    topoCmpFile=$outpath/$id.topocmp
    topoCmpLogFile=$outpath/$id.topocmp.log
    $binpath/compare_topos.py -i $msatopoSeqFile -o $topoCmpFile -log $topoCmpLogFile
 
    ### Step 6: frequency analysis
    freqAnaFile=$outpath/$id.blast.topocmp.anafreq
    echo "$SHUDATA3/wk/MPTopo/bin/histoTopocmp.py -i $topoCmpFile -bins $SHUDATA3/wk/MPTopo/test1/pidBin1.txt -pc 0 -pf 1 -o $freqAnaFile"
    $binpath/histoTopocmp.py -i $topoCmpFile -bins $SHUDATA3/wk/MPTopo/test1/pidBin1.txt -pc 0 -pf 1 -o $freqAnaFile

    ### Step 7: plotting
    echo "$SHUDATA3/wk/MPTopo/bin/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile"
    $binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile

    echo 
}
#}}}
function AnaMSATopo1() # $id#{{{
{
    local id=$1
    ### step 1: predict the topology of the query sequence using, 
    queryFastaFile=$aaPath/$id.aa
    queryTopoFile=$outpath/$id.topo
    /bin/cp -f $queryTopoFile $outpath/$id.fa
    $binpath/mySCAMPI_multi_run.pl $queryFastaFile -o $queryTopoFile
    
    ### step 2: obtain the fasta sequences for those with E-values above the threshold
    ### obtain only hits from the last round, for multiple hits from the same sequence, take the one with lower evalue
    ### sequences are obtained by fastacmd
    blastFile=$blastPath/$id.blast
    tmpFastaFile=$outpath/$id.blast.tmp.fa
    $binpath/blastm9tofasta.py $blastFile -blastdb $BLASTDB/uniprotWG100 -o $tmpFastaFile       # done

    ### Step 3: get pairwise sequence alignment of query sequence to homologues
    alignFile=$outpath/$id.blast.needle
    needle -asequence $queryFastaFile -bsequence $tmpFastaFile -gapopen 10.0 -gapextend 0.5 -outfile $alignFile

    ### Step 4: add labels of the predicted query topology to homologues
    labeledFastaFile=$outpath/$id.blast.labeled.fa
    $binpath/labeltopologyfastaseq.py -topo $queryTopoFile -aln $alignFile -fasta $tmpFastaFile -o $labeledFastaFile
 
    ### Step 5: MP topology prediction for each labeled sequence 
    ### using scampi single 
    $binpath/mySCAMPI_run.pl $labeledFastaFile --outpath $outpath
    resFile=$labeledFastaFile.res
    # problems!!!!! modhmm_scampi has bugs when using labels
    topologyFile=$outpath/$id.blast.topo
#     /scampi/compacttopo2top.py $resFile > $topologyFile
  
    ### Step 6: get topology alignment
    topologyAlignFile=$outpath/$id.topoaln
    $binpath/matchTopoPairAln.py -aln $alignFile -qtopo $queryTopoFile -ttopo $topologyFile -o $topologyAlignFile
 
    ### Step 7: pairwise comparison of topologies
    topoCmpFile=$outpath/$id.topocmp
    topoCmpLogFile=$outpath/$id.topocmp.log
    $binpath/compare_topos.py -i $topologyAlignFile  -o $topoCmpFile -log $topoCmpLogFile

    ### Step 8: obtain logodds scores for HMM prediction
    scampiOutputFile=$labeledFastaFile.xml.res
    hmmScoreFile=$outpath/$id.hmmscorelist
    $binpath/getHMMScore.py $scampiOutputFile -o $hmmScoreFile

    ### Step 9: Combine topology comparison file and hmm score file
    topoCmpHMMScoreFile=$outpath/$id.blast.topocmp_and_hmmscore
    paste $topoCmpFile $hmmScoreFile > $topoCmpHMMScoreFile
    awk '/^[^#]/{if($2!=$14) print $2, $14}' $topoCmpHMMScoreFile > $errFile
    if [ -s $errFile ]; then 
        echo "$id:  $topoCmpFile and $hmmScoreFile does not match!"  >&2
    fi

    ### Step 10: select homologs with different topologies by 
    ### a. negative DG value
    ### b. high predicted confidence 
    ### c. high homologs confidence (by evalue and pid)
    ### if pid>= 30 and logodds>= 5.8 and localCmp != OK
#     awk '/^[^#]/{if($12>=30 && $16 >= 5.8 && $6 != "OK" && $18=="yes") print}' $topoCmpHMMScoreFile | sort -k12,12rg -k16,16rg> $outpath/$id.blast.result
    $binpath/selectDIFFTopo.sh $id -binpath $binpath -dgpredpath $dgpredpath -outpath $outpath

    ### additional analysis
    ### Step add.1: frequency analysis
    freqAnaFile=$outpath/$id.blast.topocmp.anafreq
    $binpath/histoTopocmp.py -i $topoCmpFile -bins $SHUDATA3/wk/MPTopo/test1/pidBin1.txt -pc 0 -pf 1 -o $freqAnaFile

    ### Step add.2: plotting
    echo "$SHUDATA3/wk/MPTopo/bin/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile"
    $binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile
    echo 
}
#}}}
function AnaMSATopo2() # $id#{{{
{
    ###### using position specific gap opens for alignment 

    local id=$1
    ### step 1: predict the topology of the query sequence using, 
    queryFastaFile=$aaPath/$id.aa
    queryTopoFile=$outpath/$id.topo
    /bin/cp -f $queryFastaFile $outpath/$id.fa
    if [ $isPrintVerboseInfo ];then 
        echo "Step 1: $binpath/mySCAMPI_multi_run.pl $queryFastaFile -o $queryTopoFile"
    fi
    $binpath/mySCAMPI_multi_run.pl $queryFastaFile -o $queryTopoFile
    
    ### step 2: obtain the fasta sequences for those with E-values above the threshold
    ### obtain only hits from the last round, for multiple hits from the same sequence, take the one with lower evalue
    ### sequences are obtained by fastacmd
    blastFile=$blastPath/$id.blast
    homologyFastaFile=$outpath/$id.blast.fa
    if [ $isPrintVerboseInfo ];then 
        echo "Step 2: $binpath/blastm9tofasta.py $blastFile -blastdb $BLASTDB/uniprotWG100 -o $homologyFastaFile "
    fi
    $binpath/blastm9tofasta.py $blastFile -blastdb $BLASTDB/uniprotWG100 -o $homologyFastaFile

    ### Step 3: MP topology prediction for unlabelbed sequences
    ### using scampi single 
    if [ $isPrintVerboseInfo ]; then 
        echo "Step 3: $binpath/mySCAMPI_run.pl $homologyFastaFile --outpath $outpath  "
    fi
    $binpath/mySCAMPI_run.pl $homologyFastaFile --outpath $outpath
    homologyTopoFile=$homologyFastaFile.topo
    
    ### Step 4: add dg value to topology files
    if [ $isPrintVerboseInfo ];then 
        echo "Step 4: Add dg value to topology files"
    fi
    ##4.1. extract the amino acid fragment for TM regions
    homologyHelixFragFile=$outpath/$id.blast.helixfrag
    queryHelixFragFile=$outpath/$id.helixfrag
    $binpath/topo2TMfrag.py  $homologyTopoFile -f $homologyFastaFile -o $homologyHelixFragFile 
    $binpath/topo2TMfrag.py  $queryTopoFile    -f $queryFastaFile -o $queryHelixFragFile

    ##4.2. calculate DG values for all TM segments
    homologyDGFile=$outpath/$id.blast.dg
    queryDGFile=$outpath/$id.dg
    $dgpredpath/calc_dG.pl $homologyHelixFragFile  -outformat text > $homologyDGFile
    $dgpredpath/calc_dG.pl $queryHelixFragFile     -outformat text > $queryDGFile

    ##4.3. add DG values to topo files
    queryTopoWithDGScoreFile=$outpath/$id.topowithdgscore
    homologyTopoWithDGScoreFile=$outpath/$id.blast.fa.topowithdgscore
    queryDGScoreListFile=$outpath/$id.dgscorelist
    homologyDGScoreListFile=$outpath/$id.blast.fa.dgscorelist
    $binpath/topoAddDGscore.py $queryTopoFile    -dg $queryDGFile    --outpath $outpath
    $binpath/topoAddDGscore.py $homologyTopoFile -dg $homologyDGFile --outpath $outpath

    ### Step 5: Add gapopen array to Fasta File
    if [ $isPrintVerboseInfo ];then 
        echo "Step 5: Add gapopen array to fasta file"
    fi
    queryFastaWithGapOpenFile=$outpath/$id.fawithgapopen
    homologyFastaWithGapOpenFile=$outpath/$id.blast.fawithgapopen
    $binpath/fastaSeqAddGapOpenArray.py $queryFastaFile -topo $queryTopoWithDGScoreFile -o $queryFastaWithGapOpenFile
    $binpath/fastaSeqAddGapOpenArray.py $homologyFastaFile -topo $homologyTopoWithDGScoreFile -o $homologyFastaWithGapOpenFile

    ### Step 6: get pairwise sequence alignment of query sequence to homologues
    alignFile=$outpath/$id.blast.needle
    if [ $isPrintVerboseInfo ];then 
        echo "Step 6: my_needle $queryFastaWithGapOpenFile $homologyFastaWithGapOpenFile -o $alignFile  "
    fi
    my_needle $queryFastaWithGapOpenFile $homologyFastaWithGapOpenFile -o $alignFile
  
    ### Step 7: get topology alignment
    topologyAlignFile=$outpath/$id.topoaln
    if [ $isPrintVerboseInfo ];then 
        echo "Step 7: $binpath/matchTopoPairAln.py -aln $alignFile -qtopo $queryTopoFile -ttopo $homologyTopoFile -o $topologyAlignFile"
    fi
    $binpath/matchTopoPairAln.py -aln $alignFile -qtopo $queryTopoFile -ttopo $homologyTopoFile -o $topologyAlignFile
 
    ### Step 8: pairwise comparison of topologies
    if [ $isPrintVerboseInfo ];then 
        echo "Step 8: $binpath/compare_topos.py -i $topologyAlignFile  -o $topoCmpFile -log $topoCmpLogFile "
    fi
    topoCmpFile=$outpath/$id.topocmp
    topoCmpLogFile=$outpath/$id.topocmp.log
    $binpath/compare_topos.py -i $topologyAlignFile  -o $topoCmpFile -log $topoCmpLogFile

    ### Step 9: obtain logodds scores for HMM prediction
    if [ $isPrintVerboseInfo ];then 
        echo "Step 9: $binpath/getHMMScore.py $scampiOutputFile -o $hmmScoreFile "
    fi
    scampiOutputFile=$homologyFastaFile.xml.res
    hmmScoreFile=$outpath/$id.hmmscorelist
    $binpath/getHMMScore.py $scampiOutputFile -o $hmmScoreFile

    ### Step 10: Combine topology comparison file and hmm score file
    if [ $isPrintVerboseInfo ];then 
        echo "Step 10: topoCmpHMMScoreFile=$outpath/$id.blast.topocmp_and_hmmscore "
    fi
    topoCmpHMMScoreFile=$outpath/$id.blast.topocmp_and_hmmscore
    paste $topoCmpFile $hmmScoreFile > $topoCmpHMMScoreFile
    awk '/^[^#]/{if($2!=$14) print $2, $14}' $topoCmpHMMScoreFile > $errFile
    if [ -s $errFile ]; then 
        echo "$id:  $topoCmpFile and $hmmScoreFile does not match!"  >&2
    fi

    ### Step 11: select homologs with different topologies by 
    if [ $isPrintVerboseInfo ];then 
        echo "Step 11: select homologs with different topologies by DG"
    fi
    ### a. negative DG value
    ### b. high predicted confidence 
    ### c. high homologs confidence (by evalue and pid)
    ### if pid>= 30 and logodds>= 5.8 and localCmp != OK
#     awk '/^[^#]/{if($12>=30 && $16 >= 5.8 && $6 != "OK" && $18=="yes") print}' $topoCmpHMMScoreFile | sort -k12,12rg -k16,16rg> $outpath/$id.blast.result
    # add DG values to topology comparison files
    $binpath/resultAddAlnDGScore.py -topocmp $topoCmpHMMScoreFile -topoaln $topologyAlignFile  -qdg $queryDGScoreListFile -tdg $homologyDGScoreListFile -o $outpath/$id.topocmp_and_hmmscore.topoalnwithdgscore
    $binpath/selTopoAlnWithDGScore.py $outpath/$id.topocmp_and_hmmscore.topoalnwithdgscore -o $outpath/$id.sel.topoalnwithdgscore

    ### additional analysis
    ### Step add.1: frequency analysis
    if [ $isPrintVerboseInfo ];then 
        echo "Step A.1: frequency analysis"
    fi
    freqAnaFile=$outpath/$id.blast.topocmp.anafreq
    $binpath/histoTopocmp.py -i $topoCmpFile -bins $SHUDATA3/wk/MPTopo/test1/pidBin1.txt -pc 0 -pf 1 -o $freqAnaFile

    ### Step add.2: plotting
    if [ $isPrintVerboseInfo ];then 
        echo "Step A.2: plotting"
    fi
    echo "$binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile"
    $binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile
    echo 
}
#}}}
function AnaMSATopo3() # $id#{{{
{
    ###### using position specific gap opens for alignment 
    ###### the gapopen arrays are set by the sliding DG values instead of the predicted topology

    local id=$1
    ### step 1: predict the topology of the query sequence using, 
    queryFastaFile=$aaPath/$id.aa
    queryTopoFile=$outpath/$id.topo
    /bin/cp -f $queryFastaFile $outpath/$id.fa
    if [ $isPrintVerboseInfo ];then 
        echo "Step 1: $binpath/mySCAMPI_multi_run.pl $queryFastaFile -o $queryTopoFile"
    fi
    $binpath/mySCAMPI_multi_run.pl $queryFastaFile -o $queryTopoFile
    
    ### step 2: obtain the fasta sequences for those with E-values above the threshold
    ### obtain only hits from the last round, for multiple hits from the same sequence, take the one with lower evalue
    ### sequences are obtained by fastacmd
    blastFile=$blastPath/$id.blast
    homologyFastaFile=$outpath/$id.blast.fa
    if [ $isPrintVerboseInfo ];then 
        echo "Step 2: $binpath/blastm9tofasta.py $blastFile -blastdb $BLASTDB/uniprotWG100 -o $homologyFastaFile "
    fi
    $binpath/blastm9tofasta.py $blastFile -blastdb $BLASTDB/uniprotWG100 -o $homologyFastaFile

    ### Step 3: MP topology prediction for unlabelbed sequences
    ### using scampi single 
    if [ $isPrintVerboseInfo ];then 
        echo "Step 3: $binpath/mySCAMPI_run.pl $homologyFastaFile --outpath $outpath  "
    fi
     $binpath/mySCAMPI_run.pl $homologyFastaFile --outpath $outpath
    homologyTopoFile=$homologyFastaFile.topo
    
    ### Step 4: add dg value to topology files
    if [ $isPrintVerboseInfo ];then 
        echo "Step 4: Add dg value to topology files"
    fi
    ##4.1. extract the amino acid fragment for TM regions
    homologyHelixFragFile=$outpath/$id.blast.helixfrag
    queryHelixFragFile=$outpath/$id.helixfrag
    $binpath/topo2TMfrag.py  $homologyTopoFile -f $homologyFastaFile -o $homologyHelixFragFile 
    $binpath/topo2TMfrag.py  $queryTopoFile    -f $queryFastaFile -o $queryHelixFragFile

    ##4.2. calculate DG values for all TM segments
    homologyDGFile=$outpath/$id.blast.dg
    queryDGFile=$outpath/$id.dg
    $dgpredpath/calc_dG.pl $homologyHelixFragFile  -outformat text > $homologyDGFile
    $dgpredpath/calc_dG.pl $queryHelixFragFile     -outformat text > $queryDGFile

    ##4.3. add DG values to topo files
    queryTopoWithDGScoreFile=$outpath/$id.topowithdgscore
    homologyTopoWithDGScoreFile=$outpath/$id.blast.fa.topowithdgscore
    queryDGScoreListFile=$outpath/$id.dgscorelist
    homologyDGScoreListFile=$outpath/$id.blast.fa.dgscorelist
    $binpath/topoAddDGscore.py $queryTopoFile    -dg $queryDGFile    --outpath $outpath
    $binpath/topoAddDGscore.py $homologyTopoFile -dg $homologyDGFile --outpath $outpath

    #=================this step is different from  AnaMSATopo2
    ### Step 5: Add gapopen array to Fasta File
    ### Step 5.1 Get DG scan 
    homologyDGScanFile=$outpath/$id.blast.dgscan
    queryDGScanFile=$outpath/$id.dgscan
    $dgpredpath/myscanDG.pl $homologyFastaFile -lmin 21 -lmax 21 -o $homologyDGScanFile
    $dgpredpath/myscanDG.pl $queryFastaFile -lmin 21 -lmax 21 -o $queryDGScanFile

    if [ $isPrintVerboseInfo ];then 
        echo "Step 5: Add gapopen array to fasta file"
    fi
    queryFastaWithGapOpenFile=$outpath/$id.fawithgapopen
    homologyFastaWithGapOpenFile=$outpath/$id.blast.fawithgapopen
    $binpath/fastaSeqAddGapOpenArray.py $queryFastaFile -dgscan $queryDGScanFile -o $queryFastaWithGapOpenFile
    $binpath/fastaSeqAddGapOpenArray.py $homologyFastaFile -dgscan $homologyDGScanFile -o $homologyFastaWithGapOpenFile
    #=================

    ### Step 6: get pairwise sequence alignment of query sequence to homologues
    alignFile=$outpath/$id.blast.needle
    if [ $isPrintVerboseInfo ];then 
        echo "Step 6: my_needle $queryFastaWithGapOpenFile $homologyFastaWithGapOpenFile -o $alignFile  "
    fi
    my_needle $queryFastaWithGapOpenFile $homologyFastaWithGapOpenFile -o $alignFile
  
    ### Step 7: get topology alignment
    topologyAlignFile=$outpath/$id.topoaln
    if [ $isPrintVerboseInfo ];then 
        echo "Step 7: $binpath/matchTopoPairAln.py -aln $alignFile -qtopo $queryTopoFile -ttopo $homologyTopoFile -o $topologyAlignFile"
    fi
    $binpath/matchTopoPairAln.py -aln $alignFile -qtopo $queryTopoFile -ttopo $homologyTopoFile -o $topologyAlignFile
 
    ### Step 8: pairwise comparison of topologies
    if [ $isPrintVerboseInfo ];then 
        echo "Step 8: $binpath/compare_topos.py -i $topologyAlignFile  -o $topoCmpFile -log $topoCmpLogFile "
    fi
    topoCmpFile=$outpath/$id.topocmp
    topoCmpLogFile=$outpath/$id.topocmp.log
    $binpath/compare_topos.py -i $topologyAlignFile  -o $topoCmpFile -log $topoCmpLogFile

    ### Step 9: obtain logodds scores for HMM prediction
    if [ $isPrintVerboseInfo ];then 
        echo "Step 9: $binpath/getHMMScore.py $scampiOutputFile -o $hmmScoreFile "
    fi
    scampiOutputFile=$homologyFastaFile.xml.res
    hmmScoreFile=$outpath/$id.hmmscorelist
    $binpath/getHMMScore.py $scampiOutputFile -o $hmmScoreFile

    ### Step 10: Combine topology comparison file and hmm score file
    if [ $isPrintVerboseInfo ];then 
        echo "Step 10: topoCmpHMMScoreFile=$outpath/$id.blast.topocmp_and_hmmscore "
    fi
    topoCmpHMMScoreFile=$outpath/$id.blast.topocmp_and_hmmscore
    paste $topoCmpFile $hmmScoreFile > $topoCmpHMMScoreFile
    awk '/^[^#]/{if($2!=$14) print $2, $14}' $topoCmpHMMScoreFile > $errFile
    if [ -s $errFile ]; then 
        echo "$id:  $topoCmpFile and $hmmScoreFile does not match!"  >&2
    fi

    ### Step 11: select homologs with different topologies by 
    if [ $isPrintVerboseInfo ];then 
        echo "Step 11: select homologs with different topologies by DG"
    fi
    ### a. negative DG value
    ### b. high predicted confidence 
    ### c. high homologs confidence (by evalue and pid)
    ### if pid>= 30 and logodds>= 5.8 and localCmp != OK
#     awk '/^[^#]/{if($12>=30 && $16 >= 5.8 && $6 != "OK" && $18=="yes") print}' $topoCmpHMMScoreFile | sort -k12,12rg -k16,16rg> $outpath/$id.blast.result
    # add DG values to topology comparison files
    $binpath/resultAddAlnDGScore.py -topocmp $topoCmpHMMScoreFile -topoaln $topologyAlignFile  -qdg $queryDGScoreListFile -tdg $homologyDGScoreListFile -o $outpath/$id.topocmp_and_hmmscore.topoalnwithdgscore
    $binpath/selTopoAlnWithDGScore.py $outpath/$id.topocmp_and_hmmscore.topoalnwithdgscore -o $outpath/$id.sel.topoalnwithdgscore

    ### additional analysis
    ### Step add.1: frequency analysis
    if [ $isPrintVerboseInfo ];then 
        echo "Step A.1: frequency analysis"
    fi
    freqAnaFile=$outpath/$id.blast.topocmp.anafreq
    $binpath/histoTopocmp.py -i $topoCmpFile -bins $SHUDATA3/wk/MPTopo/test1/pidBin1.txt -pc 0 -pf 1 -o $freqAnaFile

    ### Step add.2: plotting
    if [ $isPrintVerboseInfo ];then 
        echo "Step A.2: plotting"
    fi
    echo "$binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile"
    $binpath/plotAnaHisto.sh -outstyle png -outpath $outpath $freqAnaFile
    echo 
}
#}}}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

errFile=/tmp/runMSATopoAna.$$.$RANDOM.err
idListFile=
aaPath=
blastPath=
outpath=
isUseOrigSeq=true #by default, use the sequence reobtained by fastacmd 

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit;;
            -l|--idlist) idListFile=$2;shift;;
            -aapath|--aapath) aaPath=$2;shift;;
            -blastpath|--blastpath) blastPath=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -notuseorig|--notuseorig) isUseOrigSeq=false;;
            -anamode|--anamode) anamode=$2;shift;;
            -v|-verbose|--verbose) isPrintVerboseInfo=1;shift;;
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
if [ ! -d "$blastPath" ]; then 
    echo "Error!  blastPath= \"$blastPath\" does not exist. Exit..."
    exit
fi
if [ ! -d "$aaPath" ]; then 
    echo "Error!  aaPath= \"$aaPath\" does not exist. Exit..."
    exit
fi

mkdir -p $outpath


for id in $(cat $idListFile); do 
    case $anamode in 
        0) AnaMSATopo0 $id;;
        1) AnaMSATopo1 $id;;
        2) AnaMSATopo2 $id;;
        3) AnaMSATopo3 $id;;
    esac
done                                                                 

rm -f $errFile

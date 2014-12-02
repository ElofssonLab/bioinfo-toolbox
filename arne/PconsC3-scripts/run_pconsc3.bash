#!/bin/bash -x


# A script to run PconsC3 and all necessary pre-calculations

# Just two variables. The location of (all) scripts and the number of cores to use (max 4 I think)
bin=/scratch/arne/PconsC3/bin
cpu=4


# Create working directory (to not fill up one directory)(

seqfile=$1
currdir=`pwd -P`

workdir=$2 
if [[ $workdir == '' ]];then 
    workdir=`echo '$dir=int(rand(100000));$dir=".PconsC3.$dir.$seqfile";if(-d $dir||-e $dir){ }else{print $dir}' | perl - $pdbid` 
    workdir=`pwd`/$workdir
fi
if [ ! -d $workdir ] 
then
    mkdir -p $workdir
    if [ $? -ne 0 ];then echo "ERROR cannot make $workdir" ; exit $? ; fi
fi
cp $seqfile $workdir

echo $seqfile $workdir
rootname=`echo $seqfile | sed -E "s/\..*//"`
cd $workdir

if [ $? -ne 0 ];then echo "ERROR cannot cd to $workdir" ; exit $? ; fi



# Make alignments
if [ ! -s $seqfile.hhE0.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE0 -e 1 $seqfile           ; fi
if [ ! -s $seqfile.hhE4.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE4 -e 1.e4 $seqfile	   ; fi
if [ ! -s $seqfile.hhE10.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE10 -e 1.e-10 $seqfile	   ; fi
if [ ! -s $seqfile.hhE40.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE40 -e 1.e-40 $seqfile	   ; fi
#
if [ ! -s $seqfile.jhE0.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE0 -e 1 $seqfile	   ; fi
if [ ! -s $seqfile.jhE4.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE4 -e 1.e4 $seqfile	   ; fi
if [ ! -s $seqfile.jhE10.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE10 -e 1.e-10 $seqfile  ; fi
if [ ! -s $seqfile.jhE40.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE40 -e 1.e-40 $seqfile  ; fi


# Obtaining statistics needed
for i in hhE0 hhE4 hhE10 hhE40 jhE0 jhE4 jhE10 jhE40
do
    if [ ! -s $seqfile.$i.stats ] ; then 
	$bin/a3mToTrimmed.py $seqfile.$i.a3m > $seqfile.$i.trimmed 
	$bin/alignmentstats.py $seqfile.$i.trimmed
    fi
done


# Running other prediction programs

# psipred can be replaced by addss.pl save some time

if [ ! -s $seqfile.ss2 ]
then   
    $bin/addss.pl $seqfile.jhE0.a3m $seqfile.jhE0.addss -a3m 
    if [ -s $seqfile.jhE0.ss2 ]
    then
	cp $seqfile.jhE0.ss2 $seqfile.ss2
    else
# This is just an extra check as addss.pl sometimes fails.
	$bin/runpsipred.py $seqfile 
	if [ -s $rootname.ss2 ]
	then
	    cp $rootname.ss2 $seqfile.ss2
	else
	    exit 1
	fi
    fi
fi



if [ ! -s $seqfile.rsa ] ; then $bin/runnetsurfp.py $seqfile ; fi 



# Running contact predictions
for i in hhE0 hhE4 hhE10 hhE40 jhE0 jhE4 jhE10 jhE40
do
    if [ ! -s $seqfile.$i.gdca ] ; then  $bin/rungdca.py $seqfile.$i.trimmed > $seqfile.$i.gneff ; fi
    if [ ! -s $seqfile.$i.plmdca20 ] ; then $bin/runplmdca20.py -c $cpu $seqfile.$i.trimmed ; fi
done

if [ ! -s $seqfile.phycmap ] 
then  
    j=`echo $seqfile | sed -E "s/\..*//g"`
    $bin/runPhyCMAP.bash  $seqfile -cpu $cpu 
    mv $j.rrunsort $seqfile.rrunsort
    mv $j.rr $seqfile.phycmap
fi


# Check that all files exist 
for i in hhE0 hhE4 hhE10 hhE40 jhE0 jhE4 jhE10 jhE40
do
    for j in gdca plmdca20 stats trimmed
    do 
	if [ ! -s $seqfile.$i.$j ] ; then echo "Missing " $seqfile.$i.$j ; exit 1  ;  fi
    done
done
for j in phycmap rsa ss2
do 
    if [ ! -s $seqfile.$j ] ; then echo "Missing " $seqfile.$j ; exit 1  ;  fi
done

# Running PconsC3


    if [ ! -s $seqfile.PconsC3.pconsc26.l6 ] ; then $bin/run_pconsc3-iterations.py $seqfile ; fi


# Cleaning up
cd $currdir

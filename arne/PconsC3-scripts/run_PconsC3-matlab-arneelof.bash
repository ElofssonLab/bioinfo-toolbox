#!/bin/bash -x

# A script to run PconsC3 and all necessary pre-calculations
seqfile=$1  # Required: Sequence file
workdir=$2  # Optional: Working directory
cpu=$3      # Optional: Number of cores (default: 1)
alnfile=$4  # Optional: Alignment file in a3m format

bin="$( cd "$( dirname "$0" )" && pwd )/bin"
echo $bin


export PATH=$PATH:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64/bin:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch/bin:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0/bin:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release/bin:/scratch/arne/PconsC3/bin/../dependencies/psipred/bin:/scratch/arne/PconsC3/bin/../dependencies/blast:/scratch/arne/PconsC3/bin/../dependencies/cd-hit-v4.5.4-2011-03-07:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch:/scratch/arne/PconsC3/bin/../dependencies/hmmer-3.1b2-linux-intel-x86_64:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release:/scratch/arne/PconsC3/bin/../dependencies/plmDCA_asymmetric_v2:/scratch/arne/PconsC3/bin/../dependencies/psipred


# Create working directory (to not fill up one directory)(
currdir=`pwd -P`
if [[ $workdir == '' ]];then 
    workdir=`echo '$dir=int(rand(100000));$dir=".PconsC3.$dir.$seqfile";if(-d $dir||-e $dir){ }else{print $dir}' | perl - $pdbid` 
    workdir=`pwd`/$workdir
fi
if [[ ! -d $workdir ]] 
then
    mkdir -p $workdir
    if [ $? -ne 0 ];then echo "ERROR cannot make $workdir" ; exit $? ; fi
fi
cp $seqfile $workdir

echo $seqfile $workdir

if [  -d `pwd`/$workdir ]; then
    workdir=`pwd`/$workdir
fi

if [[ $cpu == '' ]];then
    cpu='1'
fi


cd $workdir
seqfile=`basename $seqfile`

if [ $? -ne 0 ];then echo "ERROR cannot cd to $workdir" ; exit $? ; fi

if [[ $alnfile == '' ]];then
    i='hhE0'

    # Run HHblits
    #export HHLIB=$bin/../dependencies/hhsuite-2.0.16-linux-x86_64/lib/hh
    HHLIB=$bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch/lib/hh
    if [ ! -s $seqfile.$i.a3m ] ; then $bin/runhhblits.py -c $cpu -name $i -e 1 $seqfile ; fi
    #if [ ! -s $seqfile.hhE0.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE0 -e 1 $seqfile           ; fi
    #if [ ! -s $seqfile.hhE4.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE4 -e 1.e-4 $seqfile    ; fi
    #if [ ! -s $seqfile.hhE10.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE10 -e 1.e-10 $seqfile    ; fi
    #if [ ! -s $seqfile.hhE40.a3m ] ; then $bin/runhhblits.py -c $cpu -name hhE40 -e 1.e-40 $seqfile    ; fi

    #if [ ! -s $seqfile.jhE0.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE0 -e 1 $seqfile     ; fi
    #if [ ! -s $seqfile.jhE4.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE4 -e 1.e4 $seqfile      ; fi
    #if [ ! -s $seqfile.jhE10.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE10 -e 1.e-10 $seqfile  ; fi
    #if [ ! -s $seqfile.jhE40.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhE40 -e 1.e-40 $seqfile  ; fi
    #if [ ! -s $seqfile.jhEp3.a3m ] ; then $bin/runjackhmmer.py -c $cpu -name jhEp3 -e 1000 $seqfile  ; fi

else
    i='usr'
    cp $alnfile $seqfile.$i.a3m

fi

# Obtaining statistics needed
if [ ! -s $seqfile.$i.stats ] ; then 
$bin/a3mToTrimmed.py $seqfile.$i.a3m > $seqfile.$i.trimmed 
$bin/alignmentstats.py $seqfile.$i.trimmed
fi


# Running other prediction programs

# psipred can be replaced by addss.pl save some time

if [ ! -s $seqfile.$i.ss2 ]
then
    if [  -s $seqfile.$ss2 ]
    then
	cp $seqfile.ss2 $seqfile.$i.ss2
    else
	$bin/addss.pl $seqfile.$i.a3m $seqfile.$i.addss -a3m 
	if [ -s $seqfile.$i.ss2 ]
	then
	    cp $seqfile.$i.ss2 $seqfile.ss2
	else
            echo "addss.pl failed, running Psipred by hand"
# This is just an extra check as addss.pl sometimes fails.
	    $bin/runpsipred.py $seqfile 
	    rootname=`echo $seqfile | sed -E "s/\..*//"`
	    if [ -s $rootname.ss2 ]
	    then
		mv $rootname.ss2 $seqfile.$i.ss2
	    else
		exit 1
	    fi
	fi
    fi
fi
if [ ! -s $seqfile.ss2 ]
then
    if [  -s $seqfile.$i.ss2 ]
    then
	cp $seqfile.$i.ss2 $seqfile.ss2
    fi
fi

if [ ! -s $seqfile.rsa ] ; then $bin/runnetsurfp.py $seqfile ; fi 
if [ ! -s $seqfile.netSurf.local ] ; then cp $seqfile.rsa $seqfile.netSurf.local ; fi 


# Running contact predictions
if [ ! -s $seqfile.$i.gdca ] ; then  $bin/rungdca.py $seqfile.$i.trimmed > $seqfile.$i.gneff ; fi
if [ ! -s $seqfile.$i.plmdca20 ] ; then $bin/runplmdca20.py -c $cpu $seqfile.$i.trimmed ; fi

if [ ! -s $seqfile.phycmap ] 
then  
    #j=`echo $seqfile | sed -E "s/\..*//g"`
    $bin/runPhyCMAP.bash  $seqfile -cpu $cpu 
    mv $workdir/*.rrunsort $seqfile.rrunsort
    mv $workdir/*.rr $seqfile.phycmap
fi
if [ ! -s $seqfile.phyCMAP.results ] ; then cp $seqfile.phycmap $seqfile.phyCMAP.results; fi

# Check that all files exist 
for j in gdca plmdca20 stats trimmed
do 
	if [ ! -s $seqfile.$i.$j ] ; then echo "Missing " $seqfile.$i.$j ; exit 1  ;  fi
done
for j in phycmap rsa ss2
do 
    if [ ! -s $seqfile.$j ] ; then echo "Missing " $seqfile.$j ; exit 1  ;  fi
done


# Running PconsC3
if [ ! -s $seqfile.PconsC3.pconsc25.l5 ] ; then python $bin/runprediction-pconsc25.py $seqfile $i; fi


# Cleaning up
cd $currdir

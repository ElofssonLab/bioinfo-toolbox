#!/bin/bash -lx
#SBATCH --output=plm.%A_%a.out
#SBATCH --error=plm.%A_%a.out
#SBATCH --array=1-1000
#SBATCH -c 6
#SBATCH -t  72:00:00
#SBATCH --mem 120GB
#SBATCH -A SNIC2017-11-7

list=$1
offset=$2

pos=$(($SLURM_ARRAY_TASK_ID + $offset))

id=`tail -n+$pos $list | head -n1`

#id="pdbseq/1a3aB.fa"
#id=seq_chains/5sy1_A-5sy1_B.HH1.e-10.trimmed

name=`basename $id .trimmed`


#cd $scratch

DIR=/pfs/nobackup/home/m/mircomic/PconsC3/

# DIR=/pfs/nobackup/home/m/mircomic/PconsC3
myjulia=/pfs/nobackup/home/a/arnee/Software/julia-0d7248e2ff  # (v0.6)
#myjulia=/pfs/nobackup/home/a/arnee/Software/julia-ae26b25d43  # (v0.4)
#myjulia=/pfs/nobackup/home/a/arnee/Software/julia  # (0.6.1 compiled on keb)
#myjulia=/pfs/nobackup/home/m/mircomic/julia-ae26b25d43/
#myjulia/pfs/nobackup/home/m/mircomic/julia-ae26b25d43/

export JULIA_PKGDIR=$myjulia
export PATH=$myjulia/bin/:$PATH
#DIR=/pfs/nobackup/home/m/mircomic/PconsC3
DIR=/pfs/nobackup/home/a/arnee/git/PconsC3
sleep 2 # waiting for filesystem

dir=`dirname $id`
file=`basename $id .trimmed`
scratch=$SNIC_TMP/arnee/PLM/$dir/
mkdir -p $scratch

cd $dir


if [ ! -s $file.0.02.plm20 ]
then
    $DIR//runplm.py $file.trimmed
fi
cd ../..


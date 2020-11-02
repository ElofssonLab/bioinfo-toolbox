#!/bin/bash
#SBATCH -A SNIC2020-5-300
#SBATCH -c 5
#SBATCH -t 1:00:00
#SBATCH --array=1-3
##SBATCH -p largemem
#SBATCH --error=/home/j/juliezhu/pfs/coevolve_yeast/error/%A_%a.error
#SBATCH --output=/home/j/juliezhu/pfs/coevolve_yeast/out/%A_%a.out

#load modules
ml GCC/8.2.0-2.31.1  OpenMPI/3.1.3 ADIOS2/2.5.0-Python-3.7.2

##get the argument from the input
if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi

## param setting----------------------------------------------------------------

abspath=/home/j/juliezhu/pfs/coevolve_miipa   ##project folder path
exepath=/home/j/juliezhu/pfs/bioinfo-toolbox/rbh  ##execute script path

ecoli=/home/j/juliezhu/pfs/coevolve_yeast/pdb/seq/pdbseq.fasta   ## one fasta file of all ecoli proteins' sequences
ecolidb=/home/j/juliezhu/pfs/coevolve_yeast/pdb/seq/pdbseqDB     ## blast database of ecoli fasta file
ecoli_list=/home/j/juliezhu/pfs/coevolve_yeast/pdb/ID_pdb_trimmed.txt    ## ID list of all ecoli proteins

# result saving
rbhpath=/home/j/juliezhu/pfs/test/rbh                            ## path to save rbh files 
orthologpath=/home/j/juliezhu/pfs/test/ortholog                  ## path to save ortholog file for each ecoli protein


## other param------------------------------------------------------------------
dbpath=/home/j/juliezhu/pfs/coevolve_miipa/data_raw/blastdb_ref          ##blast database for all reference proteomes from uniprot

#idx=$SLURM_ARRAY_TASK_ID
idx=$(expr $offset + $SLURM_ARRAY_TASK_ID)
proteome_list=/home/j/juliezhu/pfs/coevolve_miipa/data_raw/refome/ID_proteomename_taxinfo.txt   ## ID list of all reference proteomes from uniprot
pname=$(sed -n ${idx}p $proteome_list | cut -d ' ' -f1)
folder=$(sed -n ${idx}p $proteome_list | cut -d ' ' -f2)


#make tmpdir
tmpdir=$(mktemp -d)
#tmpdir=/scratch/jtest
echo $tmpdir
##make dir for reciprocal best hits: forward,forward_trimmed,backward,backward_trimmed,tmp
mkdir $tmpdir/forward
mkdir $tmpdir/forward_trimmed
mkdir $tmpdir/backward
mkdir $tmpdir/backward_trimmed
mkdir $tmpdir/tmp
mkdir $tmpdir/rbh_hit
mkdir $tmpdir/rbh_hit/manual
mkdir $tmpdir/ortholog

#-------------------------------------------------------------------------------------------------------------------
### 1. reciprocal best hit process
 
## psiblast forward

cd $tmpdir/forward
psiblast -query $ecoli -db $dbpath/$pname -out $pname -num_iterations 3 -evalue 0.01 -outfmt "6 qseqid sseqid qlen slen length qcovs pident qstart qend sstart send evalue" 

## remove blank lines and lines with "Search has converged!" OBS: i use '-i' means change inplace which I could only run once.
sed -i '/^Search.*/d;/^$/d' $pname

#echo 'psiblast forward done!'

##run the forward filtering

python $exepath/blast_mainfilter.py $pname $tmpdir/forward_trimmed

#echo 'psiblast forward filtering done!'

## psiblast backward

qid=$(cat $tmpdir/forward_trimmed/$pname | cut -d '	' -f2 |cut -d '|' -f2 |sort -u)

##get the fasta file of the proteome
proteome_seq=$abspath/data_raw/refome/$folder/$pname.fasta

cd $tmpdir/backward

##iteration: blast all 'target' proteins in the file

while IFS= read -r line
do

##find the sequence of the target protein $line
pat1=".*${line}.*"
pat2="^>.*"
#echo $pat1,$pat2
sequence=$(sed "0,/${pat1}/d;/${pat2}/,\$d" $proteome_seq)

##find the info line of target protein $line
id=$(grep $line $proteome_seq)
#echo $id
#echo $sequence
##save all sequences of hits in one proteome to a fa file.
echo -e "${id}\n${sequence}" >> $tmpdir/tmp/$pname.fa
done <<< "$qid"

##psiblast back
psiblast -query $tmpdir/tmp/$pname.fa -db $ecolidb -out $pname -num_iterations 3 -evalue 0.01 -outfmt "6 qseqid sseqid qlen slen length qcovs pident qstart qend sstart send evalue" 

## remove blank lines and lines with "Search has converged!" OBS: i use '-i' means change inplace which I could only run once.
sed -i '/^Search.*/d;/^$/d' $pname

#echo 'psiblast backward done!'

## run for the backward filtering
python $exepath/blast_mainfilter.py $pname $tmpdir/backward_trimmed

#echo 'psiblast backward filtering done!'

cd $tmp
## intersection of the forward hits and backward hits; save in $tmpdir/rbh_hit
python $exepath/blast_interfwbw.py $pname $tmpdir

#echo 'intersection done!'
## filtering so that one ecoli protein has one hit in each proteome
python $exepath/parsing_blast.py $pname $tmpdir

#echo 'one hit/ecoli done!'
#----------------------------------------------------------------------------------------------

#package rbh_hit folder back 
cd $tmpdir/rbh_hit/
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf $pname.tar.zst ./*
mv $pname.tar.zst $rbhpath/

#----------------------------------------------------------------------------------------------
#cd $tmpdir/ortholog
#for m in `seq 1 5810`;do
#filename=$(sed -n ${m}p /home/j/juliezhu/pfs/coevolve_yeast/data_raw/yeast/ID_saccharomuces.txt)
#touch $filename
#done

cd $tmpdir/ortholog
filelength=$(wc -l < $ecoli_list)
for m in `seq 1 $filelength`;do
file1=$(sed -n ${m}p $ecoli_list | cut -d ' ' -f1)
file2=$(sed -n ${m}p $ecoli_list | cut -d ' ' -f2)
touch $file1
touch $file2
done


#length=5
#length=$(wc -l < $tmpdir/rbh_hit/$pname)

#for i in `seq 1 $length`;do
while IFS= read -r line;do
ecoli_name=$(echo $line | cut -d ' ' -f1 | cut -d '|' -f2)
#yeast_name=$(echo $line | cut -f1 | cut -d '|' -f2)
#yeast_row=$(sed -n ${i}p $tmpdir/rbh_hit/$pname)
#yeast_name=$(sed -n ${i}p $tmpdir/rbh_hit/$pname | cut -f1 | cut -d '|' -f2)
#echo $ecoli_name,$pname
#echo $eco_row,$eco_name
#grep $yeast_row $tmpdir/rbh_hit/$pname > $tmpfile2
#cat $tmpfile2
#sed -i "s/$/    $pname/" $tmpfile2
#echo "${yeast_row}	${pname}" >> $tmpdir/ortholog/$yeast_name
echo "${line}	${pname}" >> $tmpdir/ortholog/$ecoli_name
done < $tmpdir/rbh_hit/$pname

for i in `find . -mindepth 1 -not -empty`;do
#echo $i
ecoli_name=$(echo $i | cut -d '/' -f2)
##sed -n 1p $i
cat $i >> $orthologpath/$ecoli_name
done

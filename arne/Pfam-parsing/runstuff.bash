
#Extracting individual fasta files

bin/extract-seq.py data/escherichia_coli.fasta data/e.coli/

# Running hmmscan

 for i in data/escherichia_coli.fasta data/saccharomyces_cerevisae.fasta data/homo_sapiens.fasta  ; do j=`basename $i .fasta ` ; echo $j ;  hmmscan --cpu 4 /scratch/data/Pfam/Pfam-A.hmm $i > results/$j-pfam.out ; done &

# Running IUpred
for k in e.coli sacch homo ; do for i in data/$k/* ; do j=`basename $i .fa` ; echo $j ; iupred $i > results/iupred/$k/$j.iupred ; done ; done


# Running seg
for k in e.coli sacch homo ; do for i in data/$k/* ; do j=`basename $i .fa` ; echo $j ; seg $i > results/seg/$k/$j.seg ; done ; done

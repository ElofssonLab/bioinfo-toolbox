from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation



# Open GenBank file
handle = open('data/U00096.3.gbk', 'rU')


# For each record (mitochrodrial genome, in this case)...
for record in SeqIO.parse(handle, 'genbank') :

   # Grab the entire sequence
   seq = str(record.seq)

   # Look at all features for this record
   for feature in record.features:
      
      # If it's a CDS or rRNA...
      if feature.type == 'CDS' or feature.type == 'rRNA':

         # If it contains some attribute called 'gene' save that
         if 'gene' in feature.qualifiers:
            geneName = feature.qualifiers['gene'][0]
            
         # If it contains some attribute called 'product' save that
         elif 'product' in feature.qualifiers:
            geneName = feature.qualifiers['product'][0]
            
         # Otherwise, quit.
         else:
            print 'ERROR when parsing feature:'
            print feature.qualifiers
            quit()
         
         string=seq[feature.location.start.position:
                feature.location.end.position]
         sequence=Seq(string)
         reverse=sequence.reverse_complement()

         seqrec=SeqRecord(sequence)
         seqrec.id=geneName + '_ECOLI'
         revrec=SeqRecord(reverse)
         revrec.id=geneName + '_ECOLI'
         



         # Open output file for this gene
         geneFile = open('data/' + record.id + "/" + geneName + '_ECOLI'+ 
          '.fa', 'w')
         
         # Write FASTA header
#         geneFile.write('>')
#         geneFile.write( geneName + '_ECOLI'   + "\n")
         

         
         
         # Print this gene's seq. from complete sequence
         # Sorry about the line breaks after the next two lines
         # They should be removed, next three lines should be one.
         if feature.strand == -1:
            SeqIO.write(revrec, geneFile, "fasta")
         else:
            SeqIO.write(seqrec, geneFile, "fasta")

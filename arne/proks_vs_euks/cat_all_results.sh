cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*full* > df_full.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*domains_other* > df_domains.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*domains_shared* > df_pfam.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*linkers* > df_linkers.csv


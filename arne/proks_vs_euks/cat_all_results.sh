pushd /pfs/nobackup/home/w/wbasile/proks_euks/results/
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*full* > df_full.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*domains_other* > df_domains.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*domains_shared* > df_pfam.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*_linkers* > df_linkers.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*_N-linkers* > df_N-linkers.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*_M-linkers* > df_M-linkers.csv
cat uniprot_pfam_annotations/group_1.fasta_annotation.header uniprot_pfam_annotations/*_C-linkers* > df_C-linkers.csv
popd


#!/usr/bin/perl -w
#

use Bio::DB::Taxonomy;
use Data::Dumper;

my $outpath = "/data3/data/taxonomy/bacteria";
my $namesfile = "/data3/data/taxonomy/names.dmp";
my $nodesfile = "/data3/data/taxonomy/nodes.dmp";

my $db = Bio::DB::Taxonomy->new(
    -source => "flatfile",
    -namesfile => $namesfile,
    -nodesfile => $nodesfile,
    -directory=>"/tmp"
);

# Gram-positive bacteria correspond to Actinobacteria and Firmicutes in the NCBI Taxonomy.
# Gram-negative bacteria are all other eubacteria, except Tenericutes (including Mycoplasma), which seem to lack a type I signal peptidase and therefore do not have standard signal peptides.
#
my @gram_positive_taxname = ("Actinobacteria", "Firmicutes");
my @all_taxname = ("Actinobacteria", "Firmicutes", "eubacteria", "Tenericutes", "Mycoplasma", "Eukaryota");

my @taxidlist_actin = ();
my @taxidlist_firm = ();
my @taxidlist_eub = ();
my @taxidlist_tene = ();
my @taxidlist_myco = ();
my @taxidlist_eukar = ();

foreach my $taxname (@all_taxname){
    my $taxonid = $db->get_taxonid($taxname);
    print "$taxname => $taxonid\n";
    my $outfile = "$outpath/$taxname.taxidlist";
    my $rtvalue = GetAllTaxID_from_taxname($taxname, $outfile);
    if ($taxname eq "Actinobacteria"){
        @taxidlist_actin = @$rtvalue;
    }elsif($taxname eq "Firmicutes"){
        @taxidlist_firm = @$rtvalue;
    }elsif($taxname eq "eubacteria"){
        @taxidlist_eub = @$rtvalue;
    }elsif($taxname eq "Tenericutes"){
        @taxidlist_tene = @$rtvalue;
    }elsif($taxname eq "Mycoplasma"){
        @taxidlist_myco = @$rtvalue;
    }elsif($taxname eq "Eukaryota"){
        @taxidlist_eukar = @$rtvalue;
    }
    print "$taxname.taxidlist output"."\n";
}

print "actin = " . scalar (@taxidlist_actin)."\n";
print "firm  = " . scalar (@taxidlist_firm)."\n";
print "eub   = " . scalar (@taxidlist_eub)."\n";
print "tene  = " . scalar (@taxidlist_tene)."\n";
print "myco  = " . scalar (@taxidlist_myco)."\n";

my %hash_actin = map { $_ => 1 } @taxidlist_actin;
my %hash_firm = map { $_ => 1 } @taxidlist_firm;
my %hash_eub = map { $_ => 1 } @taxidlist_eub;
my %hash_tene = map { $_ => 1 } @taxidlist_tene;
my %hash_myco = map { $_ => 1 } @taxidlist_myco;

my %hash_gram_pos = (%hash_actin, %hash_firm);
my %hash_gram_neg = ();

foreach my $key (keys %hash_eub) {
    next if (exists $hash_actin{$key} or exists $hash_firm{$key} or exists $hash_tene{$key});
    $hash_gram_neg{$key} = 1;
}

my @gram_pos_taxidlist = keys %hash_gram_pos;
my @gram_neg_taxidlist = keys %hash_gram_neg;
print "gram_pos = ". scalar @gram_pos_taxidlist . "\n";
print "gram_neg = ". scalar @gram_neg_taxidlist . "\n";

WriteIDList(\@gram_pos_taxidlist, "$outpath/gram_pos.taxidlist");
WriteIDList(\@gram_neg_taxidlist, "$outpath/gram_neg.taxidlist");

my $gram_pn_euk_file = "$outpath/gram_pn_euk_def.txt";
open(OUT, ">", $gram_pn_euk_file) or die "Failed to write to \"$gram_pn_euk_file\"\n";
foreach my $taxid (@gram_pos_taxidlist){
    print OUT "$taxid\tgram+\n";
}
foreach my $taxid (@gram_neg_taxidlist){
    print OUT "$taxid\tgram-\n";
}
foreach my $taxid (@taxidlist_eukar){
    print OUT "$taxid\teuk\n";
}
close(OUT);


sub WriteIDList{
    my $ref_idlist = shift;
    my $outfile = shift;
    my @idlist = @$ref_idlist;
    open (OUT , ">", $outfile) or return warn "failed to write to \"$outfile\"\n";
    foreach my $taxid (@idlist){
        print OUT "$taxid\n";
    }
    close(OUT);
}

sub GetAllTaxID_from_taxname{#{{{
#usage GetAllTaxID_from_taxname(taxname, outfile)
    my $taxname = shift;
    my $outfile = shift;
    my $taxonid = $db->get_taxonid($taxname);

    my $taxon = $db->get_taxon(-taxonid => $taxonid);
    my @taxa = $db->get_all_Descendents($taxon);
    my $num = scalar(@taxa);
    my @taxidlist = ();
    open (OUT , ">", $outfile) or return;
    print OUT "$taxonid\n";
    push (@taxidlist, $taxonid);
    foreach my $tt (@taxa){
        $taxonid = $tt->ncbi_taxid();
        print OUT $taxonid."\n";
        push (@taxidlist, $taxonid);
    }
    close(OUT);
    return \@taxidlist;
}#}}}



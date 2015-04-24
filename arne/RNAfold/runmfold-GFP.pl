#!/usr/bin/perl -w

my $RUNMFOLD="/scratch/arne/RNAfold/bin/runmfold.pl";
use strict;
my $TMPDIR="/tmp/mfold".$$; 
mkdir  $TMPDIR;
#print $TMPDIR;
#chdir $TMPDIR;

my $seqdir="data/U00096.3/";
my $seqname="_ECOLI.fa";
my $csvdir="data/";
my $csvfile=$csvdir."GFP-data.csv";


open (my $CSV,$csvfile) or next;
while (my $csv=<$CSV>){
#    print $csv;
    chomp($csv);
    my @csv=split(/,/,$csv);
    my $id=lc(substr($csv[0],0,3)).uc(substr($csv[0],3,1));
    
#    print "$id\n";
#while (my $id=<>){
    my $file="data/U00096.3/".$id."_ECOLI.fa";
#    printf "Opening %s \n",$file;
    open (my $INP,$file) or 	next;
    my $seq='';
    my $inseq='';
    while (my $input=<$INP>){
#	print $input;
	chomp($input);
	if (  $input =~ /^[ACUGT][ACUGT]/){
	    $inseq.=$input;
#	    print "INPUT ",$input,"\n";
	}
    }

    my $energyA=9999;
    my $energyB=9999;
    my $energyC=9999;

    $seq='TTTAAGAAGGAGACTCGAGGAUG'.substr($inseq,3,33);
#    printf "%s %s \n",$RUNMFOLD,$seq;
    my @OUT=qx($RUNMFOLD $seq 2>/dev/null); # == 0 or die "mfold call died";
    foreach (@OUT){
	chomp;
	if (/Energy/) {
	    my @foo=split;
	    $energyA=$foo[2];
	}
    }

    $seq=substr($inseq,0,56);
#    printf "%s %s \n",$RUNMFOLD,$seq;
    @OUT=qx($RUNMFOLD $seq 2>/dev/null); # == 0 or die "mfold call died";
    foreach (@OUT){
	chomp;
	if (/Energy/) {
	    my @foo=split;
	    $energyB=$foo[2];
	}
    }

    $seq=substr($inseq,33,56);
#    printf "%s %s \n",$RUNMFOLD,$seq;
    @OUT=qx($RUNMFOLD $seq 2>/dev/null); # == 0 or die "mfold call died";
    foreach (@OUT){
	chomp;
	if (/Energy/) {
	    my @foo=split;
	    $energyC=$foo[2];
	}
    }


    printf "%s\t%s\t%s\t%f\t%f\t%f\n",$id,$seq,$energyA,$csv[8],$energyB,$energyC ;
}



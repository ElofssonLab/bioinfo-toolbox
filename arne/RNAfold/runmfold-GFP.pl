#!/usr/bin/perl -w

my $RUNMFOLD="/home/arnee/RNAfold/bin/runmfold.pl";
use strict;
my $TMPDIR="/tmp/mfold".$$; 
mkdir  $TMPDIR;
#print $TMPDIR;
#chdir $TMPDIR;


while (my $id=<>){
    chomp ($id);
    my $file="data/U00096.2/".$id."_ECOLI.fa";
#    printf "Opening %s \n",$file;
    open (my $INP,$file) or next;
    my $seq='';
    while (my $input=<$INP>){
#	print $input;
	chomp($input);
	if (  $input =~ /^[ACGT][ACGT]/){
#	    print "INPUT ",$input,"\n";
	    $seq='TTTAAGAAGGAGACTCGAGGAUG'.substr($input,4,33);
	}
    }
#    printf "%s %s \n",$RUNMFOLD,$seq;
    my @OUT=qx($RUNMFOLD $seq 2>/dev/null); # == 0 or die "mfold call died";
    foreach (@OUT){
	chomp;
	if (/Energy/) {
	    my @foo=split;
	    printf "%s\t%s\t%s\t%s\n","Energy:",$id,$seq,$foo[2] ;
	}
    }
}



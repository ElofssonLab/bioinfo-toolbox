#!/usr/bin/perl -w

my $MFOLD="/scratch/arne/RNAfold/bin/mfold-minimal.sh";
use strict;
my $TMPDIR="/tmp/mfold".$$; 
mkdir  $TMPDIR;
#print $TMPDIR;
chdir $TMPDIR;

open (my $OUT,">".$TMPDIR."/sequence") or die "cannot open OUT\n";
my $SEQ=$ARGV[0];
    printf $OUT "> Sequence\n%s\n",$SEQ;
close $OUT;

my @mfold=($MFOLD,"SEQ=sequence");
#my @OUT=system(@mfold) == 0 or die "mfold call died";
my @OUT=qx(@mfold 2>/dev/null); # == 0 or die "mfold call died";
foreach (@OUT){
    chomp;
    if (/Minimum folding energy/) {
	my @foo=split;
	printf "%s\t%s\t%s\n","Energy:",$SEQ,$foo[4] ;
    }
}

#rmdir $TMPDIR;


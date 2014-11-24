#!/usr/bin/perl -w
use strict;
my $MFOLD="/scratch/arne/RNAfold/bin/mfold-minimal.sh";
my $TMPDIR="/tmp/mfold".$$; 
mkdir  $TMPDIR;
chdir $TMPDIR;


# Definition of sequence
my %nucl = (
    A=>['A'],
    C=>['C'],
    G=>['G'],
    T=>['T'],
    Y=>['T','C'],
#    U=>['A','G'],
    U=>['U'],
    N=>['A','C','G','T']
    );


# TTTAAGAAGGAGACTCGAGGATGAGTCACTCATCCGCCCCCGAAAGGGCTACTGGAG
# TTTAAGAAGGAGACNNNNNNATGAGYCAYTCATCCGCCCCCGAAAGGGCTACTGGAG
# TTTAAGAAGGAGACNNNNNNATGUCNCAYTCATCCGCCCCCGAAAGGGCTACTGGAG
# my $seq='TTTAAGAAGGAGACNNNNNNATGAGYCAYTCATCCGCCCCCGAAAGGGCTACTGGAG';
my $seq=$ARGV[0];

my $len=length($seq);

my $i=0;
my @string='';
my @orgstring='';
my @count='';
while ($i<$len){
    $count[$i]=0;
    $orgstring[$i]=substr($seq,$i,1);
    $string[$i]=$nucl{$orgstring[$i]}[0];
    $i++;
}
$i=0;
my $j=0;

# Sequence 1

open (my $OUT,">".$TMPDIR."/sequence") or die "cannot open OUT\n";
printf $OUT "> %s\n",$j;
my $SEQ='';
foreach(@string){print $OUT $_; $SEQ.=$_;}
print $OUT "\n";
close $OUT;
my @mfold=($MFOLD,"SEQ=sequence");
#my @OUT=system(@mfold) == 0 or die "mfold call died";
my @OUT=qx(@mfold 2>/dev/null); # == 0 or die "mfold call died";
foreach (@OUT){
    chomp;
    if (/Minimum folding energy/) {
	my @foo=split;
	printf "Energy:\t%s\t%s\n",$SEQ,$foo[4] ;
    }
}
#qx('rm -rf $TMPDIR/*');


my $last=0;

my $pos=0;
while ($pos<$len){
#    my $foo=substr($seq,$pos,1);
#    print "TEST ",$pos,$count[$pos],$orgstring[$pos],"\n";
    if (exists($nucl{$orgstring[$pos]}[$count[$pos]+1]) ){
	$count[$pos]++;
	$string[$pos]=$nucl{$orgstring[$pos]}[$count[$pos]];
	for ($i=0;$i<$pos;$i++){
	    $string[$i]=$nucl{$orgstring[$i]}[0];
	    $count[$i]=0;
	}
	$j++;
	open (my $OUT,">".$TMPDIR."/sequence") or die "cannot open OUT\n";
	printf $OUT "> %s\n",$j;
	my $SEQ='';
	foreach(@string){print $OUT $_; $SEQ.=$_;}
	print $OUT "\n";
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
#	exit 0;
#	printf "> %s\n",$j;
#	foreach(@string){print $_;}
#	print "\n";
	for (my $k=0; $k<$pos; $k++){
	    $count[$k]=0;
	    $string[$k]=$nucl{substr($seq,$k,1)}[0];
	}
	$pos=0;
    }else{
	$pos++;
    }

}



#

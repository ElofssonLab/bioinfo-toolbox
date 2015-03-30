#!/usr/bin/perl -w
use strict;

my $RUNMFOLD="/home/arnee/RNAfold/bin/runmfold.pl";


my $TMPDIR="/tmp/mfold".$$; 
mkdir  $TMPDIR;
#print $TMPDIR;
#chdir $TMPDIR;

my $seqdir="data/U00096.3/";
my $seqname="_ECOLI.fa";
my $csvdir="data/";
my $csvname="_facsseq_selected.csv";

my $id=$ARGV[0];
my $file=$seqdir.$id.$seqname;
#    printf "Opening %s \n",$file;
open (my $INP,$file) or next;
my $seq='';
my %CodonsDict = ('TTT'=> 1, 'TTC'=> 2, 'TTA'=> 3, 'TTG'=> 4, 'CTT'=> 5, 
  'CTC'=> 6, 'CTA'=> 7, 'CTG'=> 8, 'ATT'=> 9, 'ATC'=> 10, 
  'ATA'=> 11, 'ATG'=> 12, 'GTT'=> 13, 'GTC'=> 14, 'GTA'=> 15, 
  'GTG'=> 16, 'TAT'=> 17, 'TAC'=> 18, 'TAA'=> 19, 'TAG'=> 20, 
  'CAT'=> 21, 'CAC'=> 22, 'CAA'=> 23, 'CAG'=> 24, 'AAT'=> 25, 
  'AAC'=> 26, 'AAA'=> 27, 'AAG'=> 28, 'GAT'=> 29, 'GAC'=> 30, 
  'GAA'=> 31, 'GAG'=> 32, 'TCT'=> 33, 'TCC'=> 34, 'TCA'=> 35, 
  'TCG'=> 36, 'CCT'=> 37, 'CCC'=> 38, 'CCA'=> 39, 'CCG'=> 40, 
  'ACT'=> 41, 'ACC'=> 42, 'ACA'=> 43, 'ACG'=> 44, 'GCT'=> 45, 
  'GCC'=> 46, 'GCA'=> 47, 'GCG'=> 48, 'TGT'=> 49, 'TGC'=> 50, 
  'TGA'=> 51, 'TGG'=> 52, 'CGT'=> 53, 'CGC'=> 54, 'CGA'=> 55, 
  'CGG'=> 56, 'AGT'=> 57, 'AGC'=> 58, 'AGA'=> 59, 'AGG'=> 60, 
  'GGT'=> 61, 'GGC'=> 62, 'GGA'=> 63, 'GGG'=> 65); 
   


while (my $input=<$INP>){
#	print $input;
    chomp($input);
    if (  $input =~ /^[ACGT][ACGT]/){
	$seq.=$input;
#	    print "INPUT ",$input,"\n";
#	    $seq='TTTAAGAAGGAGACNNNNNNATGNNNNNN'.substr($input,10); # ,27);
    }
}
my $csvfile=$csvdir.$id.$csvname;
open (my $CSV,$csvfile) or next;
while (my $csv=<$CSV>){
#	print $input;
    chomp($csv);
    my @csv=split(/,/,$csv);
    my $expr=0.0;
    my $testseq='';
    if (  $csv[0] =~ /^[ACGT][ACGT]/){
	if (  $csv[5] =~ /^[ACGT][ACGT]/){
	    $expr=$csv[7];
	}else{
	    $expr=$csv[5];
	}
	
#	    printf "%s %f\n",'TTTAAGAAGGAGAC'.$csv[0].substr($seq,6,28),$expr;
	$testseq='TTTAAGAAGGAGAC'.$csv[0].substr($seq,6,28);
	my @OUT=qx($RUNMFOLD $testseq 2>/dev/null); # == 0 or die "mfold call died";
	my $energy=9999.0;
	foreach (@OUT){
	    chomp;
	    if (/Energy/) {
		my @foo=split;
		$energy=$foo[2]
		    #printf "%s\t%s\t%s\t%s\n","Energy:",$id,$seq,$foo[2] ;
	    }
	}
	my $GC= () = $csv[0] =~ /[GC]/g;  
# codon
	my $codonA=$CodonsDict{substr($csv[0],9,3)};
	my $codonB=$CodonsDict{substr($csv[0],12,3)};
#foreach $key (sort(keys %ENV)) {
	printf "%s %s %f %f %f %d %d\n",$csv[0],$testseq,$expr,$energy,$GC/15,$codonA,$codonB;
    }
}
#    printf "%s %s \n",$RUNMFOLD,$seq;






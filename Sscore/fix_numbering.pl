#!/usr/bin/perl -w -I /proj/wallner/users/x_bjowa/git/source/perl/
#/home/denny/blandat/perl
#hej
require "bjornlib.pl";
if(scalar(@ARGV)==0) {
    print STDERR "fix_numbering.pl <model> <template> <read_seq_from_atom_in_residue (if def)>\n"; 
    exit;

}
$pdb_model=$ARGV[0];
$pdb_template=$ARGV[1];
$use_CA=1; #to get sequence
if(defined($ARGV[2])) {
    $use_CA=0;
}
$a=0;
if($use_CA) {
    ($seq_model,$resnum_model)=aa321_resnum($pdb_model);
} else {
    ($seq_model,$resnum_model)=aa321_resnumANY($pdb_model);
}
$a=1 if(defined($resnum_model));

if($use_CA) {
    ($seq_template,$resnum_template)=aa321_resnum($pdb_template);
} else {
    ($seq_template,$resnum_template)=aa321_resnumANY($pdb_template);

}
($ali_model,$ali_template)=align($seq_model,$seq_template);
@ali_model=split(//,$ali_model);
@ali_template=split(//,$ali_template);
@ali_resnum_model=(); 
my $pos=0;
my $insert_num=9000;
for(my $i=0;$i<scalar(@ali_template);$i++) {
    
    if($ali_model[$i] ne '-') { 
	if($ali_template[$i] ne '-') {
	    push(@ali_resnum_model,${$resnum_template}[$pos]);
	    #print "$ali_model[$i] $ali_template[$i] ${$resnum_template}[$pos]\n";
	} else {
	    push(@ali_resnum_model,$insert_num."X");
	    $insert_num++;
	    #print "$ali_model[$i] $ali_template[$i] 999X\n";
	}
    }
    $pos++ if($ali_template[$i] ne '-');
}
#exit;
#open(OUT,>"$pdb_model.num");
print STDERR "$ali_model\n$ali_template\n";
#exit;
my $old_resnum="whatever";
$pos=-1;
#ATOM   1226  N   GLY A 188A     20.672  19.160  17.606  1.00 26.27  
open(MODEL,$pdb_model);
while(<MODEL>) {
 #   print "A: ";
  #  print;
    if(/^ATOM/) {
#	my $atomno=substr($_, 7, 4);
#	my $atomtypeq=substr($_, 13, 3);
	my $resnum=substr($_,22,5);
	$resnum=~s/\s+//g;
	#print "$resnum $old_resnum $atomtype\n";
	if($old_resnum ne $resnum)
	{
	    $pos++;
#	    print "POS $ali_resnum_model[$pos] $_";
	}
	$old_resnum=$resnum;
	#print $pos."\n";
	if($ali_resnum_model[$pos]=~/X/) {
	    substr($_,22,5)=sprintf("%-5s",$ali_resnum_model[$pos]);
	} else {
	    substr($_,22,5)=sprintf("%5s",$ali_resnum_model[$pos]);
	}
	
    }
#    print "B: ";
    print;
}

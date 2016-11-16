#!/usr/bin/perl -w
require "bjornlib.pl";
$pdb_model=$ARGV[0];
$pdb_template=$ARGV[1];
$a=0;
($seq_model,$resnum_model)=aa321_resnum($pdb_model);
$a=1 if(defined($resnum_model));
($seq_template,$resnum_template)=aa321_resnum($pdb_template);

($ali_model,$ali_template)=align($seq_model,$seq_template);
@ali_model=split(//,$ali_model);
@ali_template=split(//,$ali_template);
@ali_resnum_model=(); 
my $pos=0;
my $insert_num=1000;
for(my $i=0;$i<scalar(@ali_template);$i++) {
    
    if($ali_model[$i] ne '-') {
	if($ali_template[$i] ne '-') {
	    push(@ali_resnum_model,${$resnum_template}[$pos]);
#	    print "$ali_model[$i] $ali_template[$i] ${$resnum_template}[$pos]\n";
	} else {
	    push(@ali_resnum_model,$insert_num."X");
	    $insert_num++;
	   # print "$ali_model[$i] $ali_template[$i] 999X\n";
	}
    }
    $pos++ if($ali_template[$i] ne '-');
}
#exit;
#open(OUT,>"$pdb_model.num");
#print STDERR "$ali_model\n$ali_template\n";
#exit;
my $old_resnum="whatever";
$pos=-1;
open(MODEL,$pdb_model);
while(<MODEL>) {
 #   print "A: ";
  #  print;
    if(/^ATOM/) {
#	my $atomno=substr($_, 7, 4);
#	my $atomtypeq=substr($_, 13, 3);
	my $resnum=substr($_,22,4);
	$resnum=~s/\s+//g;
	#print "$resnum $old_resnum $atomtype\n";
	if($old_resnum ne $resnum)
	{
	    $pos++;
	}
	$old_resnum=$resnum;
	#print $pos."\n";
	if($ali_resnum_model[$pos]=~/X/) {
	    substr($_,22,5)=sprintf("%-5s",$ali_resnum_model[$pos]);
	} else {
	    substr($_,22,4)=sprintf("%4s",$ali_resnum_model[$pos]);
	}
	
    }
#    print "B: ";
    print;
}

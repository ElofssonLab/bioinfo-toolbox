#!/usr/bin/perl -w
use strict;

print "Name\tNumContacts\tFractionMix\tFractionDiso\tScore\tScoreMix\tScoreDiso\tMeff\tDiso\tPPV\tPPVMix\tPPVDiso\n";
while(<>){
  if (/^STATs:.*average.*/){
#    print;
    chomp;
    $_=~s/\)//g;
    $_=~s/\%//g;
    $_=~s/\(//g;
    my @foo=split;
#    printf "%s\n",$_;
    printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f",$foo[1],$foo[3],$foo[4]/100,$foo[5]/100,$foo[7],$foo[8],$foo[9],$foo[11],$foo[13]/100;
      if (defined $foo[14]){
	printf "\t%f\t%f\t%f\n", $foo[15],$foo[17],$foo[19];
      }else{
	print "\tN/A\tN/A\tN/A\n";
      }
  }
}

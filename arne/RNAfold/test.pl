#!/usr/bin/perl
use strict;
use warnings;
 
my $filename = "data/U00096.2/"."zraS"."_ECOLI.fa";
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
while (my $row = <$fh>) {
  chomp $row;
  print "$row\n";
}

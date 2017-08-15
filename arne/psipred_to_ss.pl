#!/usr/bin/perl

printf ">%s\n",$ARGV[0];
while(<>){
    chomp;
    if (/^#/){next;}
    split;
    if (/^\s*\d/){
	printf "%s",$_[2];
    }
}
print "\n";
    

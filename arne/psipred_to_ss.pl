#!/usr/bin/perl

printf "> %s\n",$ARGV[0];
while(<>){
    chomp;
    @foo=split;
    if (/^\#/){next;}
    if (/^\s*\d/){
	printf "%s",$foo[2];
    }
}
print "\n";
    

#!/usr/bin/env perl

my $name;
my $top;
my $next=0;
while (<>){
    chomp;
    
    if (/^Sequence name:\s+(.*)/){
	$name=$1;
	$name=~s/\s+/\./g;
	print "> ",$name,"\n";
    }elsif (/^TOPCONS predicted topology:/){
	$next=1;
    }elsif($next){
	$top=$_;
#	print "$top\n";
	for (my $i=0;$i<length($top);$i+=60){
	    printf "%s\n",substr $top,$i,60;
	}
	$next=0;
    }
}

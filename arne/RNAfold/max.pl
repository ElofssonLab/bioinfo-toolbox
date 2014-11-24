#!/usr/bin/perl

my $max=-999990;
while(<>){
    chomp;
    my @test=split;
    my $gfp=$test[0];
    my $mg=$test[1];
    my $rna=$test[2];
    if ($max < $mg){
	printf "%s\n",$_; 
	$max=$mg;
    }
}

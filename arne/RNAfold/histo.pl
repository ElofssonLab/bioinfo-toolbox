#!/usr/bin/perl

while(<>){
    chomp;
    my @test=split;
    my $gfp=$test[1];
    my $mg=$test[2];
    my $rna=$test[16];
    my $hist=int($rna/2);
    if ($max{$hist} < $mg){$max{$hist}=$mg;}
    $sum{$hist}+=$mg;
    $num{$hist}++;
#    printf "test: %s ,  %s , %s , %s\n",$rna,$gfp,$data,$hist;
}
foreach (sort {$a <=> $b} (keys(%sum))){
    if ($_< -2){
	printf "%s\t%f\t%f\n",2*$_-1,$sum{$_}/$num{$_},$max{$_};
    }
}

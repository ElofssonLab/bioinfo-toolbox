#!/usr/bin/perl -w

my $name;
my $pdb;
while(<>){
    chomp;
    if (/.*(DP0\w+).*/){
	$name=$1;
#	print "$name\n"
    }
    if (/http\:\/\/www.pdb.org\/pdb\/explore\/explore.do\?structureId\=(\d\w\w\w)/){
	$pdb=$1;
	print "$name $pdb\n";
	system('wget http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb\&compression=NO\&structureId='.$pdb.' -O '.$name.'.pdb');
	exit 0;
    }
}





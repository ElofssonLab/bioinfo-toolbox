#!/usr/bin/perl
#
# runlga.list_list.pl script to evaluate molecules from given lists
# using LGA program and collect_PDB.pl script
#
# Author: Adam Zemla
# Current version:  (09/22/2006)
# Email: adamz@llnl.gov
#
# usage: ./runlga.list_list.pl list_1 list_2 selected_parameters
#
# Use as a selected parameters: -4 -sia -o2 -lw:3 -d:5
#                               -4 -stral
#
# in the "list_1" an additional parameter can be attached "$extra1"
# useful when running SCOP lists (for example domain range: 1mfa -aa1:1L:111L)
#
# in the "list_2" an additional parameter can be attached "$extra2"
# useful when running SCOP lists (for example domain range: 1mfa -aa2:1L:111L)
#
# lists may contain local files (for example: PDB_local/file_name)
# (by default PDB structures are located in the directory: /PDB/structures)
#
# Example of the format of the list1:
# 1o7d*            d1o7d1        class     -er1:385_B:421_B -er1:431_C:487_C
# PDB_local/file1  d1htya1       a.8.3.1   -er1:412_A:522_A
# 1wdn_A           1wdn_A_12_60  c.94.1.1  -er1:12_A:60_A
# 1ggg*            1ggg_all
# 1bve_B_5
#
# Example of the format of the list2:
# 1o7d_B
# 1gnp             1gnp_40_80    c.37.1.8  -aa2:40:80
# 1wdn_A           1wdn_A 
# 1awo___7
#

$|=1;

  $list_1_file=$ARGV[0];
  $list_2_file=$ARGV[1];
  $ARGV[0]="";
  $ARGV[1]="";

$par="@ARGV";

# subdirectory with executables
#$bin='bin';
#$bin='.';
$bin='.';

# subdirectory with results
$dirres='RESULTS/';

$res='.res';
$pdb='.pdb';
$lga='.lga';
$d='.';

@LINE=split(/\//,$list_1_file);
@R=reverse(@LINE);
$list1=$R[0];
@LINE=split(/\//,$list_2_file);
@R=reverse(@LINE);
$list2=$R[0];

$sum_res="Summary_LGA.$list1$d$list2";
system "rm -rf $sum_res ";

$i=0;
open(LIST1, "$list_1_file") || die "Can't open LIST $list1 ";
while (<LIST1>) {
  chop;
  @linefrom1=split(/\s+/,$_);
  $par1=$#linefrom1;
  $mol1_file=$linefrom1[0];
  $mol1=$mol1_file;
  $linefrom1[0]="";
  if($par1 > 0) { 
    $mol1=$linefrom1[1];
    $linefrom1[1]="";
    $linefrom1[2]="";
  }
  @LINE=split(/\//,$mol1);
  @R=reverse(@LINE);
  $mol1=$R[0];
  $mol1=~s/[\.]/\_/g;
  if($mol1 ne "") {
    $i=$i+1;
    $mollist1_file[$i]=$mol1_file;
    $mollist1[$i]=$mol1;
    $extralist1[$i]="@linefrom1";
    $extralist1[$i]=~s/-aa2:/-aa1:/g;
    $extralist1[$i]=~s/-er2:/-er1:/g;
    $extralist1[$i]=~s/-ch2:/-ch1:/g;
    $extralist1[$i]=~s/-gap2:/-gap1:/g;
  }
}
close(LIST1);
$N1=$i;
$i=0;
open(LIST2, "$list_2_file") || die "Can't open LIST $list2 ";
while (<LIST2>) {
  chop;
  @linefrom2=split(/\s+/,$_);
  $par2=$#linefrom2;
  $mol2_file=$linefrom2[0];
  $mol2=$mol2_file;
  $linefrom2[0]="";
  if($par2 > 0) { 
    $mol2=$linefrom2[1];
    $linefrom2[1]="";
    $linefrom2[2]="";
  }
  @LINE=split(/\//,$mol2);
  @R=reverse(@LINE);
  $mol2=$R[0];
  $mol2=~s/[\.]/\_/g;
  if($mol2 ne "") {
    $i=$i+1;
    $mollist2_file[$i]=$mol2_file;
    $mollist2[$i]=$mol2;
    $extralist2[$i]="@linefrom2";
    $extralist2[$i]=~s/-aa1:/-aa2:/g;
    $extralist2[$i]=~s/-er1:/-er2:/g;
    $extralist2[$i]=~s/-ch1:/-ch2:/g;
    $extralist2[$i]=~s/-gap1:/-gap2:/g;
  }
}
close(LIST2);
$N2=$i;

foreach $i2 (1 .. $N2) {  
  $mol2=$mollist2[$i2];
  $mol2_file=$mollist2_file[$i2];
  $extra2=$extralist2[$i2];
  foreach $i1 (1 .. $N1) {  
    $mol1=$mollist1[$i1];
    $mol1_file=$mollist1_file[$i1];
    $extra1=$extralist1[$i1];

    $model="$mol1$d$mol2";
    $parameters="$par $extra1 $extra2";
    $mol1_new="";
    $mol2_new="";
    while ($_="@$parameters") {
      shift;
      ($xx)=/(\S+)/;
      if($xx =~ /-mol1:/) {
        ($tmp1,$mol1_new) = split /\:/,$xx;
        $model="$mol1_new";
      }
      if($xx =~ /-mol2:/) {
        ($tmp2,$mol2_new) = split /\:/,$xx;
        $model="$mol2_new";
      }
    }
    if($mol1_new ne "" && $mol2_new ne "") {
      $model="$mol1_new$d$mol2_new";
    }

    $model=~s/\*//g;
    print "\nProcessing structures: $model $parameters \n";
    system "echo 'MOLECULE $mol1' > MOL2/$model ";
    system "$bin/collect_PDB.pl $mol1_file | grep -v '^MOLECULE ' | grep -v '^END' >> MOL2/$model ";
    system "echo 'END' >> MOL2/$model ";
    system "echo 'MOLECULE $mol2' >> MOL2/$model ";
    system "$bin/collect_PDB.pl $mol2_file | grep -v '^MOLECULE ' | grep -v '^END' >> MOL2/$model ";
    system "echo 'END' >> MOL2/$model ";
    sleep 1;
    system "$bin/lga $model $parameters > $dirres$model$res";
    sleep 1;
    chop($summaryline=`grep "^SUMMARY" $dirres$model$res | sed "s/^/$model:/"`);
    $lsumm=length($summaryline);
    if($lsumm == 0) {
#      system "rm -f $dirres$model$res ";
      print "No rotation calculated \n";
      system "echo '$model:No rotation calculated ' >> ./$sum_res ";
      system "echo 'MODEL: $model ' >> ./$sum_res\.ERRORS ";
      system "cat $dirres$model$res >> ./$sum_res\.ERRORS ";
    }
    else {
      print "$summaryline \n";
      system "echo '$summaryline ' >> ./$sum_res ";
      system "cat TMP/$model$pdb >> $dirres$model$res";
    }
    system "rm -f MOL2/$model TMP/$model$lga ";
  }
}

print "\nAll done! \n";
system "echo 'All done! (Lists: $list_1_file , $list_2_file)' >> $sum_res ";

exit;

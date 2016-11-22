#!/usr/bin/perl
#
# collect_PDB script
# script to collect PDB structures.
#
# Author: Adam Zemla 
# Email: adamz@llnl.gov
# Date: 01/26/1999
# Current version:  (09/15/2011)
#
# usage: ./collect_PDB.pl  molecule [filename]
#
# where: molecule - (input) name of the local molecule file or a code of PDB structure
#        filename - (output) name of the file where molecule will be saved (optional)
#
#        PDB code:  1abc  or  1abc_B (if chain B is requested)
#           if  1ja1_*  or  1JA1_*  then all chains are collected
#           if  1ja1_A  or  1JA1A   then chain A is collected
#           if  1JA1_   or  1ja1    then chain ' ' is collected
#
#        in the case of NMR models:
#               1bve_B_5   for PDB entry: 1bve, chain: 'B', model: 5
#               1awo___7   for PDB entry: 1awo, chain: ' ', model: 7
#
#        local files can be specified using full path:
#               DIR/molecule
#

$|=1;

$pdblocal='-';
if($#ARGV < 0) {
  print"ERROR: (Usage) collect_PDB.pl  name \n";
  exit;
}
else {
  $molecule=$ARGV[0];
  if($#ARGV > 0) { $pdblocal=$ARGV[1]; }
}

# Subdirectories where program searches for structures
$home_dir=$ENV{"HOME"};
$dirpdb="/PDB/structures";                       # default location of standard PDB entries
$dirpdbopt="PDB/structures";                     # optional location of standard PDB entries
$dirlocalpdb='PDB_local';                        # default location of local PDB structures 
$diruserspdb='PDB_users';                        # default location of users PDB structures 
#

$ok=1;
if($pdblocal eq '' || $pdblocal eq ' ' || $pdblocal eq '>'  || $pdblocal eq '>>' || $pdblocal eq '\>' || $pdblocal eq '\>\>') {
  $pdblocal='-';
}
else {
  if(-s $pdblocal > 100) { $ok=0; }
}
if($ok==1) {
  open(OUT,">$pdblocal") || {$ok=0};
  if($ok==0) {
    print "ERROR: (WRITE) bad file name: $pdblocal \n";
    exit;
  }
  else {
    $sp=' ';
    $nmr=0;
    $one='*';
    $hetatm=1;
    $l=length($molecule);
    if($molecule =~ /\// || $l>=10 || $l<4 || ($l==4 && $molecule =~ /\_/)) {
      @LINE=split(/\//,$molecule);
      @R=reverse(@LINE);
      $model=$R[0];
      $sub0="$molecule";
      $ok=1;
      if(-s $sub0 < 100) { $ok=0; } else { open(IN,$sub0) || {$ok=0}; }
    }
    else {             # checking $dirpdb and $dirlocalpdb for standard PDB files: pdbnccc.ent
      $model="$molecule";
      $pdb='pdb';
      $ent='.ent';
      $un='_';

      $two=substr($model,1,2);
      $two=~tr/A-Z/a-z/;
      $four=substr($model,0,4);
      $four=~tr/A-Z/a-z/;
      $ch5=substr($model,4,1);
      $ch7=substr($model,6,1);

      if($l>6) {
        if($ch5 eq '_' && $ch7 eq '_') {
          $one=substr($model,5,1);
        }
      }
      elsif($l==6) {
        if($ch5 eq '_' || $ch5 eq ':' ||  $ch5 eq '.') {
          $one=substr($model,5,1);
        }
      }
      else {
        if($l==5) {
          $one=substr($model,4,1);
        }
        else {
          $one=' ';
        }
      }
      if($one eq '_' || $one eq '-' || $one eq '.' || $one eq ':') {
        $one=' ';
      }
      if($l>=8) {
        $nmr=substr($model,7,$l-7);
      }
      if($nmr eq '-') {
        $nmr=0;
      }
      if($one =~ /[a-zA-Z0-9]/) { $c=$one; }
      else { $c=' '; }

      if($l<8) {
        $sub0="$dirlocalpdb/$four";            # four letter code (PDB entry) : nccc
        $sub1="$dirlocalpdb/$four$un$c";       # six letter code_chain (PDB entry): nccc_X
        $sub2="$dirlocalpdb/$pdb$four$ent";    # standard PDB files: pdbnccc.ent
        $sub3="$dirpdb/$pdb$four$ent";
        $sub4="$dirpdb/$two/$pdb$four$ent";
        $sub5="$dirpdbopt/$pdb$four$ent";
        $sub6="$dirpdbopt/$two/$pdb$four$ent";
        $sub7="$home_dir$dirpdb/$pdb$four$ent";
        $sub8="$home_dir$dirpdb/$two/$pdb$four$ent";
      }
      else {
        $sub0="$dirlocalpdb/$pdb$four$ent";    # standard PDB files: pdbnccc.ent
        $sub1="$dirpdb/$pdb$four$ent";
        $sub2="$dirpdb/$two/$pdb$four$ent";
        $sub3="$dirpdbopt/$pdb$four$ent";
        $sub4="$dirpdbopt/$two/$pdb$four$ent";
        $sub5="$home_dir$dirpdb/$pdb$four$ent";
        $sub6="$home_dir$dirpdb/$two/$pdb$four$ent";
        $sub7="$dirlocalpdb/$four";            # four letter code (PDB entry) : nccc
        $sub8="$dirlocalpdb/$four$un$c";       # six letter code_chain (PDB entry): nccc_X
      }

      $ok=1;
      if(-s $sub0 < 100) { $ok=0; } else { open(IN,$sub0) || {$ok=0}; }
      if($ok==0) {
        $ok=1;
        if(-s $sub1 < 100) { $ok=0; } else { open(IN,$sub1) || {$ok=0}; }
        if($ok==0) {
          $ok=1;
          if(-s $sub2 < 100) { $ok=0; } else { open(IN,$sub2) || {$ok=0}; }
          if($ok==0) {
            $ok=1;
            if(-s $sub3 < 100) { $ok=0; } else { open(IN,$sub3) || {$ok=0}; }
            if($ok==0) {
              $ok=1;
              if(-s $sub4 < 100) { $ok=0; } else { open(IN,$sub4) || {$ok=0}; }
              if($ok==0) {
                $ok=1;
                if(-s $sub5 < 100) { $ok=0; } else { open(IN,$sub5) || {$ok=0}; }
                if($ok==0) {
                  $ok=1;
                  if(-s $sub6 < 100) { $ok=0; } else { open(IN,$sub6) || {$ok=0}; }
                  if($ok==0) {
                    $ok=1;
                    if(-s $sub7 < 100) { $ok=0; } else { open(IN,$sub7) || {$ok=0}; }
                    if($ok==0) {
                      $ok=1;
                      if(-s $sub8 < 100) { $ok=0; } else { open(IN,$sub8) || {$ok=0}; }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    $model=~s/\*//g;

# checks for models in the current (./) , $dirlocalpdb and $diruserspdb subdirectories

    if($ok==0) {

      $sub1=$model;
      $sub2="$dirpdb/$model";
      $sub3="$home_dir$dirpdb/$model";
      $sub4="$dirlocalpdb/$model";
      $sub5="$diruserspdb/$model";
      $sub6="$diruserspdb/$model\.pdb";
      $sub7="PDB/$dirlocalpdb/$model";
      $sub8="PDB/$diruserspdb/$model";
  
      $ok=1;
      $nmr=0;
      $one='*';

      if(-s $sub1 < 100) { $ok=0; } else { open(IN,$sub1) || {$ok=0}; }
      if($ok==0) {
        $ok=1;
        if(-s $sub2 < 100) { $ok=0; } else { open(IN,$sub2) || {$ok=0}; }
        if($ok==0) {
          $ok=1;
          if(-s $sub3 < 100) { $ok=0; } else { open(IN,$sub3) || {$ok=0}; }
          if($ok==0) {
            $ok=1;
            if(-s $sub4 < 100) { $ok=0; } else { open(IN,$sub4) || {$ok=0}; }
            if($ok==0) {
              $ok=1;
              if(-s $sub5 < 100) { $ok=0; } else { open(IN,$sub5) || {$ok=0}; }
              if($ok==0) {
                $ok=1;
                if(-s $sub6 < 100) { $ok=0; } else { open(IN,$sub6) || {$ok=0}; }
                if($ok==0) {
                  $ok=1;
                  if(-s $sub7 < 100) { $ok=0; } else { open(IN,$sub7) || {$ok=0}; }
                  if($ok==0) {
                    $ok=1;
                    if(-s $sub8 < 100) { $ok=0; } else { open(IN,$sub8) || {$ok=0}; }
                    if($ok==0) {
                      print OUT "ERROR: (READ) bad model file name: $model \n";
                      close(OUT);
                      exit;
                    }
                  }
                }
              }
            }
          }
        }
      }

    }

# reads model

    if($ok==1) {
      printf OUT "MOLECULE  $model \n";
      $model_NB=0;
      $residue_name_prev="#####";
      while(<IN>) {
        chop($line=$_);
        if(/^MODEL /) {
          $model_NB=$model_NB + 1;
          ($tmp,$model_NB)=/(\S+)\s+(\S+)/;
          if($nmr == 0) {
            $nmr=$model_NB;
          }
          if($model_NB == $nmr) {
            printf OUT "$line\n";
          }
        }
        elsif(/^ATOM / || /^HETATM/) {
          $print_ok=0;
          $fixed=0;
          $a_name=substr($line,12,4);
          $a_altp=substr($line,16,1);
          $r_name=substr($line,17,3);
          $c_name=substr($line,21,1);
          $r_numb=substr($line,22,5);
          $r_insr=substr($line,26,1);
          $l=length($line);
          $left=substr($line,0,16);
          $middle=substr($line,17,9);
          $right=substr($line,27,$l-27);
          $atom="$model:$a_name:$r_name:$r_numb:$c_name";
          $atom=~s/ //g;
          $residue_name="$model:$r_name:$c_name";
          $residue_name=~s/ //g;
          $residue_numb="$model:$r_numb:$c_name";
          $residue_numb=~s/ //g;
          if($ATOMS{$atom} ne $atom) {
            if($RESNUMB{$residue_numb} eq $residue_numb && $residue_name ne $residue_name_prev) {
              if($r_insr eq ' ' && $a_altp eq ' ') {
                $r_insr='A'; $fixed=1;
              }
            }
            if($RESNUMB{$residue_numb} ne $residue_numb || $residue_name eq $residue_name_prev || $fixed == 1) {
              $line="$left$sp$middle$r_insr$right";
              if(($c_name eq $one || $one eq '*') && $model_NB == $nmr) {
                $print_ok=1;
              }
              elsif(/^HETATM/ && $hetatm == 1 && $c_name eq ' ') {
                $print_ok=1;
              }
            }
            if($print_ok==1) {
              $ATOMS{$atom}="$atom";
              printf OUT "$line\n";
              if($RESNUMB{$residue_numb} ne $residue_numb) {
                $residue_name_prev=$residue_name;
                $RESNUMB{$residue_numb}=$residue_numb;
              }
            }
          }
        }
        elsif(/^TITLE / ||
              /^LGA / ||
              /^AAMOL/ ||
              /^OBSLTE / ||
              /^DBREF / ||
              /^CONECT / ||
              /^TER / ||
              /^HET/ ||
              /^JRNL / ||
              /^FORMUL / ||
              /^MASTER / ||
              /^EXPDTA / ||
              /^AUTHOR / ||
              /^REVDAT / ||
              /^COMPND / ||
              /^SEQRES / ||
              /^SOURCE / ||
              /^KEYWDS / ||
              /^CRYST1 / ||
              /^REMARK   2[ ]+RESOLUTION/ ||
              /^REMARK  Model name: / ||
              /^REMARK  Seq: / ||
              /^REMARK  Templates: / ||
              /^REMARK  AS2TS name: / ||
              /^REMARK  AS2TS score: / ||
              /^HEADER /) {
          printf OUT "$line\n";
        }
      }
      close(IN);
      printf OUT "END \n";
    }

  }
}
close(OUT);
exit;


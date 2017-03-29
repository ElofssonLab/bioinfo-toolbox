#!/usr/bin/perl -w 
#
# ddf
#
#use lib '/home/bjornw/Research/git/source/perl/';
#use lib '/afs/pdc.kth.se/home/a/arnee/MODULES/perl5/lib/site_perl/5.6.0/i386-linux/';
#BEGIN {
#    push(@INC,'/users/bjornw/lib/bioperl/');
#    push(@INC,'/work/bjornw/lib/bioperl/');
#    push(@INC,'/Users/bjorn/Research/lib/bioperl/');
#    push(@INC,'/sw/lib/perl5/');
#}
#print join("\n",@INC);

#use Bio::Ext::Align;
use Bio::Pdb;
use Bio::AlignIO;
use Bio::SimpleAlign;
#use Bio::Tools::pSW;
use Bio::LocatableSeq;
use Bio::Seq;
use File::Temp qw/ tempfile /;

sub get_residues_with_CA
{
    my $pdbfile=shift;
    my $old_resname="undef";
    #my $old_alt_loc="undef";
    my $residue="";
    my $new_pdb="";
    #my $skip=0;
    open(PDBFILE,"$pdbfile") || die "Cannot open $pdbfile\n";
    while(<PDBFILE>)
    {
	chomp;
	if(/^ATOM/)
	{
	    my $alt_loc=substr($_,16,1);
	    my $resno=substr($_, 22, 4);
	    my $insertion_code=substr($_,26,1);
	    my $resname="$resno$insertion_code";
	    if($resname ne $old_resname && $old_resname ne "undef") #new residue
	    {
		if($residue=~/CA/)
		{
		    $new_pdb.=$residue;
		}
		$residue="";
	    }
	    if($alt_loc eq "A" || $alt_loc eq " ")
	    {
		substr($_,16,1)=" ";
		$residue.=$_."\n";
	    }
	    
	    $old_resname=$resname;
	}
    }
    if($residue=~/CA/)
    {
	$new_pdb.=$residue;
    }
#    if(not($new_pdb=~/TER/))
#    {
#	$new_pdb.="TER\n";
    #}
 #   $new_pdb.="TER\nEND\n";
    close(PDBFILE);
    return $new_pdb;
}


sub align   # Takes two strings removes all dashes and returns the alignment. 
{
    my $gap=1;
    my $ext=0.5;
    my $seq1="";
    my $seq2="";
    if(scalar(@_) > 2) {
	($seq1,$seq2,$gap,$ext)=@_;
    } else {
	($seq1,$seq2)=@_;
    }
    #print "$gap $ext\n";
    my $needle_linux="/usr/bin/needle";
    my $needle_mac="/opt/local/bin/needle";
    my $osname = $^O;
    my $input1=$seq1;
    my $input2=$seq2;
    $seq1=~s/-//g;
    $seq2=~s/-//g;
    #$seq1=remove_dashes($seq1);
    #$seq2=remove_dashes($seq2);
    $seq1=~s/\n//g;
    $seq2=~s/\n//g;
    $seq1=~s/\s+//g;
    $seq2=~s/\s+//g;
    
    my ($fh1,$file1)=tempfile("/tmp/seq.XXXXXXXX");
    my ($fh2,$file2)=tempfile("/tmp/seq.XXXXXXXX");
    my ($fh3,$file3)=tempfile("/tmp/ali.XXXXXXXX");
    close($fh3);
    print $fh1 ">seq1\n$seq1\n";
    close($fh1);
    print $fh2 ">seq2\n$seq2\n";
    close($fh2);

    if($osname eq "linux" && -e $needle_linux) {
	$needle=$needle_linux;
    }
    if($osname eq "darwin" && -e $needle_mac) {
	$needle=$needle_mac;
    }

    #print "needle -aseq $file1 -bseq $file2 -gapopen 10 -gapextend 0.5 -outfile $file3\n";
    `needle -aseq $file1 -bseq $file2 -gapopen $gap -gapextend $ext -outfile $file3 > /dev/null 2>&1`;
    #print $file3."\n";
    ($ali_return1,$ali_return2)=parse_needle_output($file3);
    `rm $file1 $file2 $file3`;
   # if(length($seq1)==1 || length($seq2)==1)
   # {
   #	my $len1=length($seq1);
   #	my $len2=length($seq2);
   #	my $len3=length($input1);
   #	my $len4=length($input2);
   #	if($len1==1)
   #	{
   #	   $dashes=length($input2)-length($input1);
   #	   $ali_return1=$input1._dashes($dashes);
   #	   $ali_return2=$input2;
   #	}
   #	if($len2==1)
   #	{
   #	   $dashes=length($input1)-length($input2);
   #	   $ali_return1=$input1;
   #	   $ali_return2=$input2._dashes($dashes);
   #	}
   #	return ($ali_return1,$ali_return2)
   # }
   #
   #
   # else
    #print "$ali_return1\n$ali_return2\n";
   
    return ($ali_return1,$ali_return2);
    
}

sub parse_needle_output
{
    my $file=shift;
    my $seq1="";
    my $seq2="";
    my $header1="";
    my $header2="";
    
    my $isFirst=1;
    open(FILE,$file);
    while(<FILE>){
	next if (/^#/);
        if(/^(.{13})(.{6}\d)\ ([\w\-]+)/){
	  #  print "header:$1, seqnro:$2, seq:$3|\n";
	    my $header = $1;
	    my $seq = $3;
	    if ($isFirst){
		$seq1.=$seq;
		$header1 = $header;
		$isFirst=0;
	    }
	    else {
		$seq2.=$seq;
		$header2 = $header;
		$isFirst=1;
	    }
        }
	
    }
    close(FILE);
    if(length($seq1) == 0) {

	print STDERR "needle 2 from the EMBOSS package needs to be installed.\n";
    }
    return($seq1,$seq2);
}

sub align_legacy   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    my $input1=$seq1;
    my $input2=$seq2;
    $seq1=~s/-//g;
    $seq2=~s/-//g;
    #$seq1=remove_dashes($seq1);
    #$seq2=remove_dashes($seq2);
    $seq1=~s/\n//g;
    $seq2=~s/\n//g;
    $seq1=~s/\s+//g;
    $seq2=~s/\s+//g;
   # if(length($seq1)==1 || length($seq2)==1)
   # {
   #	my $len1=length($seq1);
   #	my $len2=length($seq2);
   #	my $len3=length($input1);
   #	my $len4=length($input2);
   #	if($len1==1)
   #	{
   #	   $dashes=length($input2)-length($input1);
   #	   $ali_return1=$input1._dashes($dashes);
   #	   $ali_return2=$input2;
   #	}
   #	if($len2==1)
   #	{
   #	   $dashes=length($input1)-length($input2);
   #	   $ali_return1=$input1;
   #	   $ali_return2=$input2._dashes($dashes);
   #	}
   #	return ($ali_return1,$ali_return2)
   # }
   #
   #
   # else
   # {
    my ($ali_return1,$ali_return2);
    my $factory=new Bio::Tools::pSW('-matrix' => '/Users/bjorn/Research/lib/bioperl-ext/Bio/Ext/Align/blosum62.bla','-gap' => 1,'-ext' => 0);
    my $seq_obj1=Bio::Seq->new(-moltype => 'protein', -seq => $seq1, -id => "seq1");
    my $seq_obj2=Bio::Seq->new(-moltype => 'protein', -seq => $seq2, -id => "seq2");
    my $aln = $factory->pairwise_alignment($seq_obj1,$seq_obj2);
    #my $alnout = new Bio::AlignIO(-format => 'fasta',
     #                             -fh     => \*STDOUT);

    #$alnout->write_aln($aln);
#    print "align\n";
#    print $aln."\n";
#    foreach my $tmp(keys(%{$aln}))
#    {
#	print $tmp."\n";
#	print $aln{$tmp};
   # }
   # print $seq_obj2->seq();
   # print "\n";
   # print $aln,"\n";
    #my $nice_ali=$factory->align_and_show($seq_obj1,$seq_obj2,*STDOUT);
    
    ($ali_return1,$ali_return2)=fix_alignment($aln,$seq_obj1,$seq_obj2);
    #print "$ali_return1\n$ali_return2\n";
    return ($ali_return1,$ali_return2);
    
}
sub align_id   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    my $input1=$seq1;
    my $input2=$seq2;
    $seq1=~s/-//g;
    $seq2=~s/-//g;
    #$seq1=remove_dashes($seq1);
    #$seq2=remove_dashes($seq2);
    $seq1=~s/\n//g;
    $seq2=~s/\n//g;
   # if(length($seq1)==1 || length($seq2)==1)
   # {
   #	my $len1=length($seq1);
   #	my $len2=length($seq2);
   #	my $len3=length($input1);
   #	my $len4=length($input2);
   #	if($len1==1)
   #	{
   #	   $dashes=length($input2)-length($input1);
   #	   $ali_return1=$input1._dashes($dashes);
   #	   $ali_return2=$input2;
   #	}
   #	if($len2==1)
   #	{
   #	   $dashes=length($input1)-length($input2);
   #	   $ali_return1=$input1;
   #	   $ali_return2=$input2._dashes($dashes);
   #	}
   #	return ($ali_return1,$ali_return2)
   # }
   #
   #
   # else
   # {
    my ($ali_return1,$ali_return2);
    my $factory=new Bio::Tools::pSW('-matrix' => '~/lib/bioperl/Bio/Ext/Align/indent.bla','-gap' => 1,'-ext' => 0);
    my $seq_obj1=Bio::Seq->new(-moltype => 'protein', -seq => $seq1, -id => "seq1");
    my $seq_obj2=Bio::Seq->new(-moltype => 'protein', -seq => $seq2, -id => "seq2");
    my $aln = $factory->pairwise_alignment($seq_obj1,$seq_obj2);
    #$factory->align_and_show($seq_obj1,$seq_obj2,*STDOUT);
    ($ali_return1,$ali_return2)=fix_alignment($aln,$seq_obj1,$seq_obj2);
    return ($ali_return1,$ali_return2);
    
}

sub alignseq   # Takes two strings removes all dashes and returns the alignment. 
{
    my ($seq1,$seq2)=@_;
    $seq1=remove_dashes($seq1);
    $seq2=remove_dashes($seq2);
    my ($ali_return1,$ali_return2);
    my $factory=new Bio::Tools::pSW('-matrix' => '/Users/bjorn/Research/lib/bioperl-ext/Bio/Ext/Align/blosum62.bla','-gap' => 12,'-ext' => 4);
    my $seq_obj1=Bio::Seq->new(-moltype => 'protein', -seq => $seq1, -id => "seq1");
    my $seq_obj2=Bio::Seq->new(-moltype => 'protein', -seq => $seq2, -id => "seq2");
    my $aln = $factory->pairwise_alignment($seq_obj1,$seq_obj2);
    #$factory->align_and_show($seq_obj1,$seq_obj2,*STDOUT);
    ($ali_return1,$ali_return2)=fix_alignment($aln,$seq_obj1,$seq_obj2);
    #$seq_obj1->DESTROY;
    #$seq_obj2->DESTROY;
    #$factory->DESTROY;
    #delete($seq_obj1);
    #delete($seq_obj2);
    #delete($factory);
    return ($ali_return1,$ali_return2);
}


sub remove_ends
{
    #$ali2 has to be the pdbfile.

    my ($ali1,$ali2)=@_;
    my $start=0;
    my $end=0;
    my $find_start=1;
    my $find_end=0;
    my @ali1 = split('',$ali1);
    my @ali2 = split('',$ali2);
    my $ali_return1=$ali1;
    my $ali_return2=$ali2;
    my @seq_2_ali=();
    for(my $i=0;$i<scalar @ali1;$i++)
    {
	push(@seq_2_ali,$i) if($ali1[$i] ne "-");
	if($ali1[$i] ne "-" && $ali2[$i] ne "-" && $find_start)
	{
	    $start=$i;
	    $find_start=0;
	    $find_end=1;
	}
	if($ali1[$i] ne "-" && $ali2[$i] ne "-")  #alignad target
	{
	    $end=$i;
	}
    }
    
    #print "\n$start $end\n";
    my $number=0;
    my $len=length($ali1);
        
    ###  REMOVE THE RESIDUES WHICH ARE NOT ALIGNED *****
    $start++;
    $end++;
    substr($ali_return1,0,$start-1)=_dashes($start-1);
    substr($ali_return1,$end,$len-$end)=_dashes($len-$end);
    my @ali_return1=split('',$ali_return1);
    my @ali_return2=split('',$ali_return2);
    my $ali_return1_1="";
    my $ali_return2_2="";
    #Remove all "double dashes"
    for(my $i=0; $i < length($ali_return1);$i++)
    {
	if(not($ali_return1[$i] eq "-" && $ali_return2[$i] eq "-"))
	{
	    $ali_return1_1.=$ali_return1[$i];
	    $ali_return2_2.=$ali_return2[$i];
	}
	#if($ali2[$i] ne "-")
	#{
	#    $start_pdb++ if($i<$start);
	#    $end_pdb++ if($i<$end);
	#}
    }
    
    #print "$ali_return1_1\n\n$ali_return2_2\n\n";
    return($ali_return1_1,$ali_return2_2);
}

sub merge_ali   # $A aligned with $B1, $B2 aligned with $C, returns $A aligned with $C.
{               # $B1 and $B2 have exactly the same residues.
    my($A,$B1,$B2,$C) = @_;
    my $A_return="";
    my $C_return="";
    my @A=split('',$A);
    my @B1=split('',$B1);
    my @B2=split('',$B2);
    my @C=split('',$C);
    my $len1 = scalar @A; #length of first aligment pair 
    my $len2 = scalar @B2; #length of second aligment pair

    my $j=0;   # $i indexes in the first alignment $j in the second
    my $i=0;
    my $temp="";;
    my $temp2="";
    while($i<$len1 || $j<$len2) 
    {
	#print "$B1[$i] $B2[$j]\n";
	if(defined($B1[$i]) && $B1[$i] eq "-")
	{
	    $temp2.="-";
	    $temp.=$A[$i];
	    $i++;
	}
	elsif($B2[$j] eq "-")
	{
	    $temp2.=$C[$j];
	    $temp.="-";
	    $j++;
	}
	elsif($B1[$i] eq $B2[$j])
	{
	    #print "$i $j $A[$i] $C[$j]\n";
	    $temp2.=$C[$j];
	    $temp.=$A[$i];
	    $i++;
	    $j++;
	}
	
    }
    my $temp3="";
    my $temp4="";
    my @temp=split('',$temp);
    my @temp2=split('',$temp2);
    #print $temp."\n\n".$temp2."\n\n";
    # Remove all positions which has gaps in both sequences.
    for(my $i=0;$i < scalar @temp2;$i++)
    {
	if(not ($temp[$i] eq "-" && $temp2[$i] eq "-"))
	{
	    $temp3.=$temp[$i];
	    $temp4.=$temp2[$i];
	}

    }

    #print $temp3."\n\n".$temp4."\n\n";
    return($temp3,$temp4);
}

sub merge_ali_do_not_remove_gaps   # $A aligned with $B1, $B2 aligned with $C, returns $A aligned with $C.
{               # $B1 and $B2 MUST have exactly the same residues!!!!!!
    my($A,$B1,$B2,$C) = @_;
    my $A_return="";
    my $C_return="";
    my @A=split('',$A);
    my @B1=split('',$B1);
    my @B2=split('',$B2);
    my @C=split('',$C);
    my $len1 = scalar @A; #length of first aligment pair 
    my $len2 = scalar @B2; #length of second aligment pair

    my $j=0;   # $i indexes in the first alignment $j in the second
    my $i=0;
    my $temp="";;
    my $temp2="";
    my $counter=0;
#    exit;

    #print "A:  $A\nB1: $B1\nB2: $B2\nC:  $C\n";
    while($i<$len1 || $j<$len2) 
    {
	#print "$B1[$i] $B2[$j]\n";
	#exit;
#	print "$B1[$i] $B2[$j] $i $len1 $j $len2\n";
#	print $counter."\n";
	if(defined($B1[$i]) && defined($B2[$j]) && $B1[$i] eq "-" && $B2[$j] eq "-")
	{
	    $temp2.=$C[$j];
	    $temp.=$A[$i];
	    $i++;
	    $j++;
	}
	elsif(defined($B1[$i]) && $B1[$i] eq "-")
	{
	    $temp2.="-";
	    $temp.=$A[$i];
	    $i++;
	}
	elsif($B2[$j] eq "-")
	{
	    $temp2.=$C[$j];
	    $temp.="-";
	    $j++;
	}
	elsif($B1[$i] eq $B2[$j])
	{
	    #print "$i $j $A[$i] $C[$j]\n";
	    $temp2.=$C[$j];
	    $temp.=$A[$i];
	    $i++;
	    $j++;
	}
	
    }
    #my $temp3="";
    #my $temp4="";
    #my @temp=split('',$temp);
    #my @temp2=split('',$temp2);
    ##print $temp."\n\n".$temp2."\n\n";
    ## Remove all positions which has gaps in both sequences.
    #for(my $i=0;$i < scalar @temp2;$i++)
    #{
    #	if(not ($temp[$i] eq "-" && $temp2[$i] eq "-"))
    #	{
    #	    #print "hej\n";
    #	    $temp3.=$temp[$i];
    #	    $temp4.=$temp2[$i];
    #	}
    #
    #}

    #print $temp3."\n\n".$temp4."\n\n";
    return($temp,$temp2);
}

sub remove_dashes
{
    my $temp=shift;
    my @list=split(/\-+/,$temp);
    my $new_str=join('',@list);
    return $new_str;

}
sub fix_alignment    # $seq1 and $seq2 must be in the same order as in the $aln otherwise it gets wrong.
{
    my($aln,$seq1,$seq2)=@_;
    my $ali_seq1="";
    my $ali_seq2="";
    # Parse the alignment
    my ($seq,@start, @end,@ali_seqs);
    my $i=0;

    foreach $seq ($aln->each_seq())
    {
	#print $seq->seq()."\n";
	$start[$i]=$seq->start();
	$end[$i]=$seq->end();
	$ali_seqs[$i]=$seq->seq();
	#print $seq."\n";
	#print $i."\n";
	#print $ali_seqs[$i]."\n";

	#print "$start[$i] $end[$i]\n";
	$i++;
    }
    #print $start[0]."\n".$start[1]."\n";
# Reformat alignment so that it contain all resides and -.
    
  #  print "$ali_seqs[0]\n$ali_seqs[1]\n";
  #  exit;
# Fix the begining
    if($start[0]!=1)
    {
	$ali_seq1.=$seq1->subseq(1,$start[0]-1);
	$ali_seq2.=_dashes($start[0]-1);
    }
    if($start[1]!=1)
    {
	$ali_seq1.=_dashes($start[1]-1);
	$ali_seq2.=$seq2->subseq(1,$start[1]-1);
    }
    
# Add the alignment
    
    $ali_seq1.=$ali_seqs[0];
    $ali_seq2.=$ali_seqs[1];
    
# Fix alignment end;
 #   print "$end[0] ".$seq1->length()."\n";
    if($end[0]<$seq1->length())
    {
	my $len=$seq1->length();
	$ali_seq1.=$seq1->subseq($end[0]+1,$len);
	$ali_seq2.=_dashes($len-$end[0]);
    }
    #print "$end[1] ".$seq2->length()."\n";
    if($end[1]<$seq2->length())
    {
	my $len=$seq2->length();
	$ali_seq1.=_dashes($len-$end[1]);
	$ali_seq2.=$seq2->subseq($end[1]+1,$len);
    }
    return ($ali_seq1,$ali_seq2);
}

sub _stars
{

my $number=shift;
    my $str="";

    for(my $i=0;$i<$number;$i++)
    {
	$str.="*";
    }
    return $str;
}
sub _dashes
{
    my $number=shift;
    my $str="";

    for(my $i=0;$i<$number;$i++)
    {
	$str.="-";
    }
    return $str;

}
sub merge_multi_ali
{
   # print "In merger\n";
    #my (@target_ali_list,@template_ali_list)=@_;
    my $temp=scalar @_;                                    #The first half of the vector is target
    my @target_ali_list=@_[0..$temp/2-1];                  #the second is template.
    my @template_ali_list=@_[$temp/2..$temp-1];
	
    #print @target_ali_list;
    #print "\n";
    #print @template_ali_list;
    #print "\n";
    my $target_return="";
    my @template_return=();
    my $res_count=0;
    my @res_found=();
    my $number_of_seq=scalar @target_ali_list;
    my @index_vec=();
    my $total_ali_length=0;
    my @temp_str=();

    # Sort the target sequences so the one with the most number of dashes in the begin is first.
    # Remember to put them in the same original order again, when done!!!
    for(my $i=0;$i<$number_of_seq;$i++)
    {
	#print "$template_ali_list[$i]\n\n$target_ali_list[$i]\n\n";
	$temp_str[$i]=count_begin_dashes($target_ali_list[$i])." ".$i." ".$template_ali_list[$i]." ".$target_ali_list[$i]." ";
    }
    @temp_str=sort numerically_str @temp_str;
    my @old_index=();
    my $trash="";
    #Initialize the vectors. And split the sequences.
    for(my $i=0;$i<$number_of_seq;$i++)
    {
	($trash,$old_index[$i],$template_ali_list[$i],$target_ali_list[$i])=split(/\s+/,$temp_str[$i]);
	push(@index_vec,0);
	push(@res_found,0);
	#push(@target_return,"");
	push(@template_return,"");
	$total_ali_length+=length $target_ali_list[$i];#  if(length $target_ali_list[$i]>$longest_ali);
	#print $target_ali_list[$i]."\n\n".$template_ali_list[$i]."\n\n";
	$target_ali_list[$i]=[split('',$target_ali_list[$i])];
	$template_ali_list[$i]=[split('',$template_ali_list[$i])];

    } 
    #exit;
    if($number_of_seq>0)
    {
	while(sum(@index_vec)<$total_ali_length)
	{
	    #print $total_ali_length."\n";
	    #print sum(@index_vec)."\n";
	    #print sum(@res_found)."\n";
	    for(my $i=0;$i<$number_of_seq;$i++)
	    {
		my $summa=sum(@res_found);
		my $summa2=sum(@index_vec);
		#print "@res_found $number_of_seq $summa $summa2 $total_ali_length $index_vec[$i] $target_ali_list[$i][$index_vec[$i]] $template_ali_list[$i][$index_vec[$i]]\n";
		if(sum(@res_found) == $number_of_seq)  #all @index_vec points a target reside align or not does not matter print 'em all..
		{
		    if(defined($target_ali_list[0][$index_vec[0]]))
		    {
			$target_return.=$target_ali_list[0][$index_vec[0]];
		    }
		    else
		    {
			$target_return.="-";
		    }

		    #print $target_return."\n";
		    for(my $j=0;$j<$number_of_seq;$j++)
		    {
			if(defined($template_ali_list[$j][$index_vec[$j]]))
			{
			    $template_return[$j].=$template_ali_list[$j][$index_vec[$j]];
			}
			else
			{
			    $template_return[$j].="-";
			}
			$res_found[$j]=0;
			$index_vec[$j]++;
		    }
		    
		}
		elsif(defined($target_ali_list[$i][$index_vec[$i]]) &&
		      defined($template_ali_list[$i][$index_vec[$i]]) &&
		      $target_ali_list[$i][$index_vec[$i]] eq "-" &&
		      $template_ali_list[$i][$index_vec[$i]] ne "-") # unaligned template res. Add that res
		                                                  # to return and "-" to all others which are not unaligned.
		{
		    for(my $j=0;$j<$number_of_seq;$j++)
		    {
			if(defined($template_ali_list[$j][$index_vec[$j]]) && $i==$j)
			{
			    $template_return[$j].=$template_ali_list[$j][$index_vec[$j]];
			}
			elsif(defined($template_ali_list[$j][$index_vec[$j]]) && 
			      defined($target_ali_list[$j][$index_vec[$j]]) &&
			      $target_ali_list[$j][$index_vec[$j]] eq "-") #if unaligned res
			{
			    $template_return[$j].=$template_ali_list[$j][$index_vec[$j]];
			    $index_vec[$j]++;
			}
			else
			{
			    $template_return[$j].="-";
			}
		    }
		    $target_return.="-";
		    $index_vec[$i]++;
		}
		elsif(defined($target_ali_list[$i][$index_vec[$i]]) &&  $target_ali_list[$i][$index_vec[$i]] ne "-") #target residue
		{
		    $res_found[$i]=1;
		}
	    }
	}
    }
    # Put the sequences back in order they came in.
    my @template_return2=();
    for(my $i=0;$i<scalar @old_index;$i++)
    {
	$template_return2[$old_index[$i]]=$template_return[$i];
	
    }

    my @temp_target=split(//,$target_return);
    my @temp_template=();
    my $i=0;
    foreach my $line (@template_return2)
    {
	 my @temp=split(//,$line);
	 $temp_template[$i]=[split('',$template_return2[$i])];
	 $i++;
    }

    
    my $start=0;
    my $end=0;
    my $find_start=1;
    my $find_end=0;
    my @seq_2_ali=();
    for(my $i=0;$i<scalar @temp_target;$i++)
    {
	push(@seq_2_ali,$i) if($temp_target[$i] ne "-");
	 if($temp_target[$i] ne "-" && $find_start)
	 {
	     my $align=0;
	     for(my $j=0;$j<scalar @temp_template;$j++)
	     {
		 $align=1 if($temp_template[$j][$i] ne "-");

	     }
	     if($align)
	     {
		 $start=$i;
		 $find_start=0;
		 $find_end=1;

	     }
	 }
	 if($temp_target[$i] ne "-")
	 {
	     my $align=0;
	     for(my $j=0;$j<scalar @temp_template;$j++)
	     {
		 $align=1 if($temp_template[$j][$i] ne "-");

	     }
	     $end=$i if($align);
	 }
    }

    
    #foreach my $str(@template_return2)
    #{
#	print $str."\n";
#    }
    #print "\n$target_return\n";
    #print "$start $end\n";
    my $len=length($target_return);
    $start++;
    $end++;


    #OBS WILL MODEL ALL RESIDUES, ADDED DURING CASP SEASON 2004-06-20 /BW
 
    #substr($target_return,0,$start-1)=_dashes($start-1);
    #substr($target_return,$end,$len-$end)=_dashes($len-$end);

    
    #remove all gaps aligned to gaps in all sequences which might have been introduced by the dashes above.
    @temp_target=split('',$target_return);
    

    for(my $i=0;$i<$number_of_seq;$i++)
    {
	$temp_template[$i]=[split('',$template_return2[$i])];

    } 
    for(my $i=0;$i<scalar @temp_target;$i++)
    {
	if($temp_target[$i] eq '-')
	{
	    my $dash_counter=1;
	    #store all residues in one column
	    for(my $j=0;$j<$number_of_seq;$j++)
	    {
		$dash_counter++ if($temp_template[$j][$i] eq '-');
	    }
	    
	    if($dash_counter == $number_of_seq+1)   #We have all gaps
	    {
		splice(@temp_target,$i,1);
		for(my $j=0;$j<$number_of_seq;$j++)
		{
		    splice(@{$temp_template[$j]},$i,1);
		}
		$i--;
		#print $dash_counter." ".$number_of_seq."\n";

	    }
	}
    }
    
    $target_return=join('',@temp_target);
    for(my $i=0;$i<$number_of_seq;$i++)
    {
	$template_return2[$i]=join('',@{$temp_template[$i]});

    } 

    #print "\n$target_return\n";

#*#    print $target_return."\n";
#*#    for(my $i=0;$i<scalar @template_return2;$i++)
#*#    {
#*#	   print "=====\n";
#*#	   print $template_return2[$i]."\n";
#*#	   
#*#    }
    return ($target_return,@template_return2);

}

sub merge_multi_ali2
{
   # print "In merger\n";
    #my (@target_ali_list,@template_ali_list)=@_;
    my $temp=scalar @_;                                    #The first half of the vector is target
    my @target_ali_list=@_[0..$temp/2-1];                  #the second is template.
    my @template_ali_list=@_[$temp/2..$temp-1];
	
    #print @target_ali_list;
    #print "\n";
    #print @template_ali_list;
    #print "\n";
    my $target_return="";
    my @template_return=();
    my $res_count=0;
    my @res_found=();
    my $number_of_seq=scalar @target_ali_list;
    my @index_vec=();
    my $total_ali_length=0;
    my @temp_str=();

    # Sort the target sequences so the one with the most number of dashes in the begin is first.
    # Remember to put them in the same original order again, when done!!!
    for(my $i=0;$i<$number_of_seq;$i++)
    {
	#print "$template_ali_list[$i]\n\n$target_ali_list[$i]\n\n";
	$temp_str[$i]=count_begin_dashes($target_ali_list[$i])." ".$i." ".$template_ali_list[$i]." ".$target_ali_list[$i]." ";
    }
    @temp_str=sort numerically_str @temp_str;
    my @old_index=();
    my $trash="";
    #Initialize the vectors. And split the sequences.
    for(my $i=0;$i<$number_of_seq;$i++)
    {
	($trash,$old_index[$i],$template_ali_list[$i],$target_ali_list[$i])=split(/\s+/,$temp_str[$i]);
	push(@index_vec,0);
	push(@res_found,0);
	#push(@target_return,"");
	push(@template_return,"");
	$total_ali_length+=length $target_ali_list[$i];#  if(length $target_ali_list[$i]>$longest_ali);
	print $target_ali_list[$i]."\n\n".$template_ali_list[$i]."\n\n";
	$target_ali_list[$i]=[split('',$target_ali_list[$i])];
	$template_ali_list[$i]=[split('',$template_ali_list[$i])];

    } 
    #exit;
    if($number_of_seq>0)
    {
	while(sum(@index_vec)<$total_ali_length)
	{
	    #print $total_ali_length."\n";
	    #print sum(@index_vec)."\n";
	    #print sum(@res_found)."\n";
	    for(my $i=0;$i<$number_of_seq;$i++)
	    {
		my $summa=sum(@res_found);
		my $summa2=sum(@index_vec);
		#print "@res_found $number_of_seq $summa $summa2 $total_ali_length $index_vec[$i] $target_ali_list[$i][$index_vec[$i]] $template_ali_list[$i][$index_vec[$i]]\n";
		if(sum(@res_found) == $number_of_seq)  #all @index_vec points a target reside align or not does not matter print 'em all..
		{
		    if(defined($target_ali_list[0][$index_vec[0]]))
		    {
			$target_return.=$target_ali_list[0][$index_vec[0]];
		    }
		    else
		    {
			$target_return.="-";
		    }

		#    print $target_return."\n";
		    for(my $j=0;$j<$number_of_seq;$j++)
		    {
			if(defined($template_ali_list[$j][$index_vec[$j]]))
			{
			    $template_return[$j].=$template_ali_list[$j][$index_vec[$j]];
			}
			else
			{
			    $template_return[$j].="-";
			}
			$res_found[$j]=0;
			$index_vec[$j]++;
		    }
		    
		}
		elsif(defined($target_ali_list[$i][$index_vec[$i]]) &&
		      defined($template_ali_list[$i][$index_vec[$i]]) &&
		      $target_ali_list[$i][$index_vec[$i]] eq "-" &&
		      $template_ali_list[$i][$index_vec[$i]] ne "-") # unaligned template res. Add that res
		                                                  # to return and "-" to all others which are not unaligned.
		{
		    for(my $j=0;$j<$number_of_seq;$j++)
		    {
			if(defined($template_ali_list[$j][$index_vec[$j]]) && $i==$j)
			{
			    $template_return[$j].=$template_ali_list[$j][$index_vec[$j]];
			}
			elsif(defined($template_ali_list[$j][$index_vec[$j]]) && 
			      defined($target_ali_list[$j][$index_vec[$j]]) &&
			      $target_ali_list[$j][$index_vec[$j]] eq "-") #if unaligned res
			{
			    $template_return[$j].=$template_ali_list[$j][$index_vec[$j]];
			    $index_vec[$j]++;
			}
			else
			{
			    $template_return[$j].="-";
			}
		    }
		    $target_return.="-";
		    $index_vec[$i]++;
		}
		elsif(defined($target_ali_list[$i][$index_vec[$i]]) &&  $target_ali_list[$i][$index_vec[$i]] ne "-") #target residue
		{
		    $res_found[$i]=1;
		}
	    }
	}
    }
    # Put the sequences back in order they came in.
    my @template_return2=();
    for(my $i=0;$i<scalar @old_index;$i++)
    {
	$template_return2[$old_index[$i]]=$template_return[$i];
	
    }

    my @temp_target=split(//,$target_return);
    my @temp_template=();
    my $i=0;
    foreach my $line (@template_return2)
    {
	 my @temp=split(//,$line);
	 $temp_template[$i]=[split('',$template_return2[$i])];
	 $i++;
    }

    
    my $start=0;
    my $end=0;
    my $find_start=1;
    my $find_end=0;
    my @seq_2_ali=();
    for(my $i=0;$i<scalar @temp_target;$i++)
    {
	push(@seq_2_ali,$i) if($temp_target[$i] ne "-");
	 if($temp_target[$i] ne "-" && $find_start)
	 {
	     my $align=0;
	     for(my $j=0;$j<scalar @temp_template;$j++)
	     {
		 $align=1 if($temp_template[$j][$i] ne "-");

	     }
	     if($align)
	     {
		 $start=$i;
		 $find_start=0;
		 $find_end=1;

	     }
	 }
	 if($temp_target[$i] ne "-")
	 {
	     my $align=0;
	     for(my $j=0;$j<scalar @temp_template;$j++)
	     {
		 $align=1 if($temp_template[$j][$i] ne "-");

	     }
	     $end=$i if($align);
	 }
    }

    
    #foreach my $str(@template_return2)
    #{
#	print $str."\n";
#    }
    #print "\n$target_return\n";
    #print "$start $end\n";
    my $len=length($target_return);
    $start++;
    $end++;
    substr($target_return,0,$start-1)=_dashes($start-1);
    substr($target_return,$end,$len-$end)=_dashes($len-$end);

    
    #remove all gaps aligned to gaps in all sequences which might have been introduced by the dashes above.
    @temp_target=split('',$target_return);
    

    for(my $i=0;$i<$number_of_seq;$i++)
    {
	$temp_template[$i]=[split('',$template_return2[$i])];

    } 
    for(my $i=0;$i<scalar @temp_target;$i++)
    {
	if($temp_target[$i] eq '-')
	{
	    my $dash_counter=1;
	    #store all residues in one column
	    for(my $j=0;$j<$number_of_seq;$j++)
	    {
		$dash_counter++ if($temp_template[$j][$i] eq '-');
	    }
	    
	    if($dash_counter == $number_of_seq+1)   #We have all gaps
	    {
		splice(@temp_target,$i,1);
		for(my $j=0;$j<$number_of_seq;$j++)
		{
		    splice(@{$temp_template[$j]},$i,1);
		}
		$i--;
		#print $dash_counter." ".$number_of_seq."\n";

	    }
	}
    }
    
    $target_return=join('',@temp_target);
    for(my $i=0;$i<$number_of_seq;$i++)
    {
	$template_return2[$i]=join('',@{$temp_template[$i]});

    } 

    #print "\n$target_return\n";

#*#    print $target_return."\n";
#*#    for(my $i=0;$i<scalar @template_return2;$i++)
#*#    {
#*#	   print "=====\n";
#*#	   print $template_return2[$i]."\n";
#*#	   
#*#    }
    return ($target_return,@template_return2);

}

sub MSA
{
    my $temp=scalar @_;                                    #The first half of the vector is target
    my @target_ali_list=@_[0..$temp/2-1];                  #the second is template.
    my @template_ali_list=@_[$temp/2..$temp-1];
    my ($target_ali_reformat,@templates)=merge_multi_ali(@target_ali_list,@template_ali_list);
    
    my @target_seq=split(//,$target_ali_reformat);
    my @temp_vec=();
    my @template_return=();
    my $target_return="";

    

    #splitting
    foreach my $seq(@templates)
    {
	my @temp=split(//,$seq);
	#print $seq."\n";
	push(@temp_vec,[@temp]);
	push(@template_return,"");
    }
    #print scalar @template_return;
    #exit;
    for(my $i=0;$i<scalar @target_seq;$i++)
    {
	 #if($target_seq[$i] ne '-' || $target_seq[$i] eq '-')
	 #{
	     $target_return.=$target_seq[$i];
	     for($j=0;$j<scalar @template_return;$j++)
	     {
		 $template_return[$j].=$temp_vec[$j][$i];
		# print "$i $j\n";
	     }

	 #}

     }
    return($target_return,@template_return);
    
}

sub histogram
{
    my ($ranges,$x,$frac)=@_;
    my @ranges=@{$ranges};
    my @x=@{$x};
    my @hist=();
    my @middle_points=();
   # print scalar @ranges," ",scalar @x," ",scalar @{$x[0]},"\n";;
    if(scalar @ranges == 1)
    {
	print STDERR "\@ranges only contain 1 element!\n";
	return 0;
    } 
   # print @{$x[0]},"\n";

   
    #print scalar @x."\n";
 #   exit;
    for(my $i=0;$i<scalar @ranges-1;$i++)
    {
	#my @temp=();
#	push(@temp,0);
	#if(scalar @x>0)
#	{
	for(my $j=0;$j<scalar @x;$j++)
	{
	#    push(@temp,0);
	
#	}
	#push(@hist,\@temp);
	    ${$hist[$i]}[$j]=0
	}
#	
    }
    
    for(my $i=0;$i<scalar @ranges - 1;$i++)
    {
	push(@middle_points,($ranges[$i]+$ranges[$i+1])/2);
	#print scalar @x,"\n===\n";
	for(my $j=0;$j<scalar @x;$j++)
	{
	   # print scalar @{$x[$j]},"\n";
	    #if(defined(@{$x[$j]}))
	    {
		for(my $k=0;$k<scalar @{$x[$j]};$k++)
		{
		    ${$hist[$i]}[$j]++ if(${$x[$j]}[$k]>=$ranges[$i] && ${$x[$j]}[$k]<$ranges[$i+1]);
		    
		    #print "${$x[$j]}[$k] $ranges[$i+1] i=$i j=$j k=$k\n";# if(not(defined($ranges[$i+1])));
		}
	    }
	    #print "\n" if($j==0);
	}
	#exit;
    }
    if(defined($frac) && $frac==1)
    {
	$len=scalar @x;
	#print $len,"\n";
	for(my $i=0;$i<scalar @ranges-1;$i++)
	{
	    my @temp=();
	    #printf("%6.2f ",$middle_points[$i]);
	    for(my $j=0;$j<scalar @x;$j++)
	    {
		#printf("before %4d ",${$hist[$i]}[$j]);
		${$hist[$i]}[$j]/=scalar(@{$x[$j]});
	#	printf("after %4d ",${$hist[$i]}[$j]/$len);
	    }
	 #   print "\n";
	}
	# $len=scalar @x;
    }


    return([@middle_points],[@hist]);

}
sub corrcoef
{
    my ($x,$y)=@_;
    my @x=@{$x};
    my @y=@{$y};
    my $xm=mean(@x);
    my $ym=mean(@y);
    if(scalar @y == 0)
    {
	warn "You need to use bracket (\[\]) for both array arguments!\n";
    }
    my $cov=0;
    my $corr=0;
    my $elements=scalar(@x);
    for(my $i=0;$i<scalar @x;$i++)
    {
#	print "$x[$i] $y[$i]\n";
	$cov+=($x[$i]-$xm)*($y[$i]-$ym)
    }
    $cov=$cov/($elements-1);
    $std_x=std(@x);
    $std_y=std(@y);
    $corr=0;
    if($std_x*$std_y != 0)
    {
	$corr=$cov/($std_x*$std_y);
    }
    return $corr;
}

sub rmsd
{
    my ($x,$y)=@_;
    my @x=@{$x};
    my @y=@{$y};
    my $sum=0;


    
    #my $xm=mean(@x);
    #my $ym=mean(@y);
    if(scalar @y == 0)
    {
	warn "You need to use bracket (\[\]) for both array arguments!\n";
    }
    #my $cov=0;
    #my $corr=0;
    my $elements=scalar(@x);
    for(my $i=0;$i<scalar @x;$i++)
    {
	$square=($x[$i]-$y[$i])*($x[$i]-$y[$i]);
#	print "($x[$i]-$y[$i])*($x[$i]-$y[$i]) $square\n";
	$sum+=($x[$i]-$y[$i])*($x[$i]-$y[$i]);
#	print "$x[$i] $y[$i]\n";
#	$cov+=($x[$i]-$xm)*($y[$i]-$ym)
    }
#    $cov=$cov/($elements-1);
#    $std_x=std(@x);
#    $std_y=std(@y);
#    $corr=0;
#    if($std_x*$std_y != 0)
#    {
#	$corr=$cov/($std_x*$std_y);
 #   }
    return sqrt($sum/$elements);
}
sub rmsd_scaled
{
    my ($x,$y)=@_;
    my @x=@{$x};
    my @y=@{$y};
    my $sum=0;


    
    #my $xm=mean(@x);
    #my $ym=mean(@y);
    if(scalar @y == 0)
    {
	warn "You need to use bracket (\[\]) for both array arguments!\n";
    }
    #my $cov=0;
    #my $corr=0;
    my $elements=scalar(@x);
    for(my $i=0;$i<scalar @x;$i++)
    {
#	$square=($x[$i]-$y[$i])*($x[$i]-$y[$i]);
#	print "($x[$i]-$y[$i])*($x[$i]-$y[$i]) $square\n";
	$x_s=1/($x[$i]*$x[$i]/5+1);
	$y_s=1/($y[$i]*$y[$i]/5+1);

	$sum+=($x_s-$y_s)*($x_s-$y_s); #1/(abs(($x[$i]-$y[$i]))+1); #*($x[$i]-$y[$i]);
#	print "$x[$i] $y[$i]\n";
#	$cov+=($x[$i]-$xm)*($y[$i]-$ym)
    }
#    $cov=$cov/($elements-1);
#    $std_x=std(@x);
#    $std_y=std(@y);
#    $corr=0;
#    if($std_x*$std_y != 0)
#    {
#	$corr=$cov/($std_x*$std_y);
 #   }
 #   my $mean=$sum/$elements;
#    $mean=1/$mean-1;

    return sqrt($sum/$elements);
}
sub error
{
    my @data=@_;
    my $mean=mean(@data);
    my $n=scalar @data;
    if($n!=1 && $n != 0)
    {
	my $sum=0;
	foreach my $term(@data)
	{
	    $sum+=($term-$mean)*($term-$mean);
	}
	return 1.960*(sqrt(1/($n*($n-1))*$sum));
    }
    else
    {
	return 0;
    }
}


sub std
{
    my @data=@_;
    my $mean=mean(@data);
    my $n=scalar @data;
    if($n!=1)
    {
	my $sum=0;
	foreach my $term(@data)
	{
	    $sum+=($term-$mean)*($term-$mean);
	}
	return sqrt(1/($n-1)*$sum);
    }
    else
    {
	return 0;
    }
}



sub mean
{
    my @data=@_;
    my $sum=0;
    my $number_of_elements=scalar @data;
    foreach my $term(@data)
    {
	$sum+=$term;
    }
    if($number_of_elements==0)
    {
	return 0;
    }
    else
    {
	return $sum/$number_of_elements;
    }
}

sub median
{
    my @data=@_;
    my @data_s=sort numerically (@data);
    my $middle=int(scalar(@data)/2);
    return($data_s[$middle]);
}


sub numerically
{
    $b<=>$a;
}

sub numerically_inc
{
    $a<=>$b
}

sub sum
{
    my @vec=@_;
    my $sum=0;
    foreach my $term(@vec)
    {
	#print $term,"\n";
	$sum+=$term;
    }
    return $sum;
}

sub max
{
    my @vec=@_;
    my $max=$vec[0];
    foreach my $term(@vec)
    {
	if($term>$max)
	{
	    $max=$term;
	}
#	print $term,"\n";
    }
    return $max;
}
sub min
{
    my @vec=@_;
    my $min=$vec[0];
    foreach my $term(@vec)
    {
	if($term<$min)
	{
	    $min=$term;
	}
	#print $term,"\n";
    }
    return $min;
}

sub shift_score
{
    
    my @vec=@_;
    my $sum=0;
    foreach my $term(@vec)
    {
	#print $term,"\n";
	if($term ne "?")
	{
	    $sum+=(1/(1+abs($term)));
	}
    }
    return $sum;

}

sub count_begin_dashes
{
    my $str=shift;
    my @temp=split('',$str);
    my $counter=0;
    for(my $i=0;$i<scalar @temp;$i++)
    {
	if($temp[$i] eq '-')
	{
	    $counter++;
	}
	else
	{
	    last;
	}
	
    }
    return $counter;

}
sub numerically_str
{
    my @temp_a=split(/\s+/,$a);
    my @temp_b=split(/\s+/,$b);
    $temp_b[0]<=>$temp_a[0];
}
#####my $special_restraints=generate_top_code($ss_pred,$target_ali_new,$target_seq);

sub generate_top_code
{
    #Assumes that the $ali_seq is a subset of $target_seq;
    
    my ($pred,$ali_seq,$target_seq)=@_;
    my($ali_seq2,$target_seq2)=align($ali_seq,$target_seq);

    #print "$ali_seq2\n\n$target_seq2\n\n$pred\n\n";
    #$ali_seq2.='-';

    my @ali_seq=split('',$ali_seq2);
    my ($begin,$end,$start,$stop)=(0,0,0,0);
    my $get_stop=0;
    #my @target_seq=split('',$target_seq2);
    for(my $i=0;$i<scalar @ali_seq;$i++)
    {
	if($ali_seq[$i] ne "-" && not($get_stop))
	{
	    $start=$i+1;
	    $get_stop=1;
	}
	if($get_stop && $ali_seq[$i] eq "-")
	{
	    $stop=$i; #+1-1
	    last;
	}

    }
    $stop=length($target_seq) if($stop==0);
    #my $len=$stop-$start+1;
    #print "$start $end $len\n";
    
    my $alpha="";
    my $strand="";
    my @pred=split(//,$pred);
    
    my $last_pred="";
    my $j=1;
    for(my $i=$start-1;$i<=$stop-1;$i++,$j++)
    {
	
	if($last_pred ne $pred[$i] && $i>$start-1)
	{
	    if($last_pred eq "H")  #print alpha restraints
	    {
		$end=$j-1;
		$alpha.="  MAKE_RESTRAINTS RESTRAINT_TYPE = 'ALPHA', RESIDUE_IDS = '$begin' '$end'\n";
	    }
	    elsif($last_pred eq "E") #print strand restraints
	    {
		$end=$j-1;
		$strand.="  MAKE_RESTRAINTS RESTRAINT_TYPE = 'STRAND',RESIDUE_IDS = '$begin' '$end'\n";
		    
	    }
	}
	if($pred[$i] eq "H" && $last_pred ne "H")  #new helix
	{
	    $begin=$j;
	}
	elsif($pred[$i] eq "E" && $last_pred ne "E") #new strand
	{
	    $begin=$j;
	}
	$last_pred=$pred[$i];	
    }
    $j--;
    if($last_pred eq "H" && $begin>$end)  #print alpha restraints, second check means that a new begin must have been added since last time..
    {
	$end=$j-1;
	
	$alpha.="  MAKE_RESTRAINTS RESTRAINT_TYPE = 'ALPHA', RESIDUE_IDS = '$begin' '$end'\n";
    }
    elsif($last_pred eq "E" && $begin>$end) #print strand restraints
    {
	$end=$j-1;
	$strand.="  MAKE_RESTRAINTS RESTRAINT_TYPE = 'STRAND',RESIDUE_IDS = '$begin' '$end'\n";
    }
    
    my $return_str="SUBROUTINE ROUTINE = 'special_restraints'\n  SET ADD_RESTRAINTS = on\n$alpha$strand";
    $return_str.="  RETURN\nEND_SUBROUTINE\n";
            
    #print $return_str;
    return $return_str;

}
sub get_start_number
{
    my ($seq1,$seq2)=@_;
    my ($target_seq,$model_seq)=align($seq1,$seq2);
    #align_loc($seq1,$seq2);
    my @res=split(//,$model_seq);
    my $start=1;
    for(my $i=0;$i<scalar @res;$i++)
    {
	if($res[$i] ne '-')
	{
	    $start=$i+1;
	    last;
	}
    }
    #print $start."\n";
    return($start);
}
sub number_pdb #Takes a pdbfile and a number to start with, and starts renumbering the file from there
{
    my($file,$number)=@_;
    $in_number=$number;
    $number--;
    my $oldresno="";
    my $old_chain=" ";
    my $temp="";
    my $new_file="";
    my $atomcount=1;
    open(PDB,"$file");
    while(<PDB>)
    {
	my $line=$_;
	if($line=~/^ATOM/)
	{
	    my $chain=substr($line,21,1);
	    my $resno=substr($line, 22, 4);
	    if($oldresno ne $resno)
	    {
		$number++;
		if($old_chain ne $chain) { 
		    $number=$in_number;
	#	    $new_file.="TER\n";
		}
		$temp=sprintf("%4d",$number);
	    }
	    substr($line,6,5)=sprintf("%5d",$atomcount);
	    $atomcount++;
	    substr($line, 22, 4)=$temp;
	    substr($line, 26, 1)=" ";
	    $oldresno=$resno;
	    $old_chain=$chain;
	}
	$new_file.=$line;
    }
    return $new_file;
}
sub read_in_psipred
{
    my $file=shift;
    #my $target_seq=shift;
    my $seq="";
    my $ss_pred="";
    open(FILE,"$file");# or die "Cannot open $file.\n";
    while(<FILE>)
    {
	if(/^Pred:\s[HCE]+/)
	{
	    my @temp=split(/\s+/);
	    $ss_pred.=$temp[1] if(scalar @temp >=2);
	}
	elsif(/^  AA:/)
	{
	    my @temp=split(/\s+/);
	    $seq.=$temp[2] if(defined($temp[2]));
	}
    }
    close(FILE);
    #print length($seq)," ",length($ss_pred),"\n";
    #print $file,"\n";
    if(length($seq)==0 && length($ss_pred)==0)
    {
	#print $file,"\n";
	open(FILE,"$file");# or die "Cannot open $file.\n";
	while(<FILE>)
	{
	    if(/\s+\d+/)
	    {
		chomp;
		@temp=split(/\s+/);
		$ss_pred.=$temp[3];
		$seq.=$temp[2];
	    }
	}
    }
    close(FILE);
    #print "$seq\n$ss_pred\n";
    #exit;
    #generate_top_code($ss_pred,10,110);
    return ($seq,$ss_pred);

}

sub get_ss
{
    my ($seq,$seq_pred,$pred)=@_;
    my $return_ss="";
    #my ($seq_pred,$pred)=read_in_psipred($psifile);
    my ($ali_real,$ali_pred)=align($seq,$seq_pred);
   # print "$ali_real\n$ali_pred\n";#$pred\n";
    #print "------------\n";
    my @ali_real=split(//,$ali_real);
    my @ali_pred=split(//,$ali_pred);
    my @pred=split(//,$pred);
    my $ss_pred_ali="";
    my $i=0;
    foreach my $residue(@ali_pred)
    {
	if($residue eq '-')
	{
	    $ss_pred_ali.="-";
	}
	else
	{
	    $ss_pred_ali.=$pred[$i];
	    $i++;
	}
    }
    my @ss_pred_ali=split(//,$ss_pred_ali);
    for(my $i=0;$i<=$#ali_real;$i++)
    {
	if($ali_real[$i] ne '-')
	{
	    $return_ss.=$ss_pred_ali[$i];
	    
	}
    }
    #print $return_ss."\n";
    return $return_ss;
}

sub get_rsa
{
    my ($seq,$seq_pred,$pred)=@_;
    my @pred=@{$pred};
    my @return_rsa=();
    #my ($seq_pred,$pred)=read_in_psipred($psifile);
    my ($ali_real,$ali_pred)=align($seq,$seq_pred);
    #print "$ali_real\n$ali_pred\n";#$pred\n";
    #print "------------\n";
    my @ali_real=split(//,$ali_real);
    my @ali_pred=split(//,$ali_pred);
  #  my @pred=split(//,$pred);
    my @ss_pred_ali=();
    my $i=0;
    foreach my $residue(@ali_pred)
    {
    	if($residue eq '-')
    	{
    	    push(@ss_pred_ali,"-");
    	}
    	else
    	{
	    push(@ss_pred_ali,$pred[$i]);
    	    $i++;
    	}
    }
    #my @ss_pred_ali=split(//,$ss_pred_ali);
    for(my $i=0;$i<=$#ali_real;$i++)
    {
	if($ali_real[$i] ne '-')
	{
	    push(@return_rsa,$ss_pred_ali[$i]);
	}
    }
    #print $return_ss."\n";
    return @return_rsa;
}

sub read_in_stride
{
    my $file=shift;
    my $seq="";
    my $ss="";
    open(FILE,"$file") || die "Cannot open $file. (bjornlib)\n";
    while(<FILE>)
    {
	chomp;
	if(/^SEQ/)
	{
	    #print substr($_,10,50)."\n";
	    $seq.=substr($_,10,50);
	}
	if(/^STR/)
	{
	    #print substr($_,10,50)."\n";
	    $ss.=substr($_,10,50);
	}
	last if(/^LOC/);
    }
    #Remove white spaces at end.
    $seq=~s/ //g;
    $ss=substr($ss,0,length($seq));
    my @temp=split(//,$ss);
    my $return_ss="";
    foreach my $ss_a(@temp)
    {
	#print $ss_a."\n";
	if($ss_a eq 'H' || $ss_a eq 'G' || $ss_a eq 'I')
	{
	    $return_ss.="H";
	}
	elsif($ss_a eq "E")
	{
	    $return_ss.="E";
	}
	else
	{
	    $return_ss.="C";
	}
    }
    #print "$seq\n$return_ss\n$ss\n";
    return($seq,$return_ss);
}

sub read_in_stride_ext
{
    my $file=shift;
    my $seq="";
    my $ss="";
    open(FILE,"$file") || die "Cannot open $file. (bjornlib)\n";
    while(<FILE>)
    {
	chomp;
	if(/^SEQ/)
	{
	    #print substr($_,10,50)."\n";
	    $seq.=substr($_,10,50);
	}
	if(/^STR/)
	{
	    #print substr($_,10,50)."\n";
	    $ss.=substr($_,10,50);
	}
	last if(/^LOC/);
    }
    #Remove white spaces at end.
    $seq=~s/ //g;
    $ss=substr($ss,0,length($seq));
    $ss=~s/ /C/g;
    $ss=~s/T/C/g;
    return($seq,$ss);
}


sub get_pdb
{
    my $pdbcode=shift;
    $pdbcode=lc($pdbcode);
    #print STDERR "PDBCODE to get_pdb $pdbcode\n";
    my $PDBURL="ftp://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/";
    my $OBSOLETE_PDBURL="ftp://ftp.rcsb.org/pub/pdb/data/structures/obsolete/pdb/";
    my $MODEL_PDBURL="ftp://ftp.rcsb.org/pub/pdb/data/structures/models/current/pdb/";
    my $PDBDIR="/afs/pdc.kth.se/projects/sbc/mirror/mirrors/pdb/";
    my $PDBDIR_OBS="/afs/pdc.kth.se/projects/sbc/mirror/mirrors/pdb/obsolete/";
    my $PDBDIR_MODELS="/afs/pdc.kth.se/projects/sbc/mirror/mirrors/pdb/models/current/pdb/";
    my $PDBDIR_MODELS_OBS="/afs/pdc.kth.se/projects/sbc/mirror/mirrors/pdb/models/obsolete/pdb/";
    my $SCOPDIR="/afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.63/pdb/";
    my $SCOPDIR2="/afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.57/pdb/";
    if(length($pdbcode)<6)
    {
	my $file="pdb$pdbcode.ent";
	my $subdir=lc(substr($pdbcode,1,2));
	my $return_file="";;
	if(-e $file)
	{
	    return $file;
	}
	#print "Trying afs directories\n";
	if(-e "$PDBDIR$subdir/$file")
	{
	    $return_file="$PDBDIR$subdir/$file";
	    #`cp $PDBDIR$subdir/$file .`;
	}
	if(!-e "$return_file" && -e "$PDBDIR_OBS$subdir/$file")
	{
	    
	    $return_file="$PDBDIR_OBS$subdir/$file";
	    #`cp $PDBDIR_OBS$subdir/$file .`;
	}
	if(!-e "$return_file" && -e "$PDBDIR_MODELS$subdir/$file")
	{
	    $return_file="$PDBDIR_MODELS$subdir/$file";
	    #`cp $PDBDIR_MODELS$subdir/$file .`;
	}
	if(!-e "$return_file" && -e "$PDBDIR_MODELS_OBS$subdir/$file")
	{
	    $return_file="$PDBDIR_MODELS_OBS$subdir/$file";
	    #`cp $PDBDIR_MODELS$subdir/$file .`;
	}
	if(!-e "$return_file")
	{
	    my $pwd=`pwd`;
	    print $pwd."\n";
	    print "Trying ordinary URL $file.Z\n";
	    my $count=0;
	    while(!-e "$file.Z" && $count<3)
	    {
		 print "Try $count: wget $PDBURL$file.Z\n";
		`wget $PDBURL$file.Z`;
		$count++;
	    }
	    if(-e "$file.Z")
	    {
		`gunzip -f $file.Z`;
	    }
	}
	if(!-e "$file" && !-e "$return_file")
	{
	    #print "Trying theoretical models URL\n";
	    my $pdbfile_compressed=$file.".Z";
	    my $pdbfile_dir="$subdir/$pdbfile_compressed";
	    $count=0;
	    while(not(-e "$file.Z") && $count<3)
	    {
		   print "Try $count: $MODEL_PDBURL$pdbfile_dir\n";
		`wget $MODEL_PDBURL$pdbfile_dir`;
		$count++;
	    }
	    if(-e "$file.Z")
	    {
		`gunzip -f $file.Z`;
	    }
	}
	
	if(!-e "$file" && !-e "$return_file")
	{
	    #print "Trying obsolete URL\n";
	    #my $subdir=substr($pdb_code,1,2);
	    $subdir.="/";
	    my $pdbfile_compressed=$file.".Z";
	    my $pdbfile_obsolete=$subdir.$pdbfile_compressed;
	    $count=0;
	    while(not(-e "$file.Z") && $count<3)
	    {
		#    print "Try $count: $OBSOLETE_PDBURL$pdbfile_obsolete\n";
		`wget $OBSOLETE_PDBURL$pdbfile_obsolete`;
		$count++;
	    }
	    if(-e "$file.Z")
	    {
		`gunzip -f $file.Z`;
	    }
	}
	if(-e "$file")
	{
	    return $file;
	}
	elsif(-e "$return_file")
	{
	    return $return_file;
	}
	else
	{
	    return 0;
	}
    }
    else #SCOP FILE
    {
	
	if(length($pdbcode)==7)
	{
	    
	    #print "egrep $pdbcode /afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.63/data/dir.cla.scop.txt_1.63|head -n 1|awk '{print \$4}'\n";
	    my $str=`egrep $pdbcode /afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.63/data/dir.cla.scop.txt_1.63|head -n 1|awk '{print \$4}'`;
	    chomp($str);
	    #print "STR $str\n";
	    $pdbcode.=".$str";
		
	}
	my $subdir=substr($pdbcode,8,1);
	#print $pdbcode."\n";
	my $file=$SCOPDIR.$subdir."/$pdbcode.pdb";
	my $file2=$SCOPDIR2.$subdir."/$pdbcode.pdb";
	#print $file."\n";
	if(-e $file)
	{
	    return $file
	}
	elsif(-e $file2)
	{
	    return $file2;
	}
	else
	{
	    return 0;
	}


    }


}



sub seq_ok
{
    my $seq1=shift;
    my $seq2=shift;
    $seq1=~s/-//g;
    $seq1=~s/A//g;
    $seq1=~s/C//g;
    $seq1=~s/G//g;
    $seq1=~s/T//g;

    $seq2=~s/-//g;
    $seq2=~s/A//g;
    $seq2=~s/C//g;
    $seq2=~s/G//g;
    $seq2=~s/T//g;
    
    my $bool=1;
    $bool=0 if(length($seq1)<2 || length($seq1)<2);
    return $bool;
}


sub id
{
    my ($seq1,$seq2)=@_;
   # my ($ali1,$ali2)=alignseq($seq1,$seq2);
    my ($ali1,$ali2)=align($seq1,$seq2,10,0.5);
    
    my $shortest_length=length($seq1);
    $shortest_length=length($seq2) if($shortest_length>length($seq2));
    print "$ali1\n$ali2\n";
    my @ali1=split(//,$ali1);
    my @ali2=split(//,$ali2);
    my $id_residues=0;
    for(my $i=0;$i<scalar @ali1;$i++)
    {
	
	$id_residues++ if($ali1[$i] eq $ali2[$i]);
    }
    #print "$ali1\n$ali2\n";
    #print "$id_residues\n";
    return $id_residues/$shortest_length;
}
sub seq_id_for_aligned_residues
{
    my ($ali1,$ali2)=@_;
    #my ($ali1,$ali2)=alignseq($seq1,$seq2);
    my $len=0;
    #my $shortest_length=length($seq1);
    #$shortest_length=length($seq2) if($shortest_length>length($seq2));
    #print "$ali1\n$ali2\n";
    my @ali1=split(//,$ali1);
    my @ali2=split(//,$ali2);
    my $id_residues=0;
    for(my $i=0;$i<scalar @ali1;$i++)
    {
	if($ali1[$i] ne "-" && $ali2[$i] ne "-")
	{
	    $id_residues++ if($ali1[$i] eq $ali2[$i]);
	    $len++;
	}
    }
    #print "$ali1\n$ali2\n";
    #print "$id_residues\n";
    return $id_residues/$len;
}

sub id_ali
{
    my ($ali1,$ali2)=@_;
    #my ($ali1,$ali2)=alignseq($seq1,$seq2);
    #my $shortest_length=length($seq1);
    #$shortest_length=length($seq2) if($shortest_length>length($seq2));
    #print "$ali1\n$ali2\n";
    my @ali1=split(//,$ali1);
    my @ali2=split(//,$ali2);
    my $id_residues=0;
    my $len1=0;
    my $len2=0;
    for(my $i=0;$i<scalar @ali1;$i++)
    {
	
	$id_residues++ if($ali1[$i] eq $ali2[$i] && $ali1[$i] ne '-');
	$len1++ if($ali1[$i] ne '-');
	$len2++ if($ali2[$i] ne '-');
    }
    my $shortest_length=$len1;
    $shortest_length=$len2 if($len2<$len1);
#print "$ali1\n$ali2\n";
   # print "$id_residues $len1 $len2\n";
    return sprintf("%5.4f",$id_residues/$shortest_length);
}
sub sspred_model
{
   

    my ($model_seq,$target_seq,$ss_pred)=@_;
    my ($model_ali,$target_ali)=align($model_seq,$target_seq);
    
    my @ss_pred=split(//,$ss_pred);
    my @model_ali=split(//,$model_ali);
    my @target_ali=split(//,$target_ali);
    my $sspred_model="";
    my $j=0;
    for(my $i=0;$i<scalar @target_ali;$i++)
    {
	if($target_ali[$i] ne "-" && $model_ali[$i] ne "-")
	{
	    $sspred_model.=$ss_pred[$j];
	}
	$j++ if($target_ali[$i] ne "-");

    }
    return $sspred_model;
}

sub readfile
{
    my $file=shift;
    open(FILE,$file);
    my @out=<FILE>;
    my $out=join('',@out);
    return $out;

}

sub aa321CA
{
    my $file=shift;
    my $seq="";
    my $old_resnum="whatever";
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 12, 4);
	    my $resnum=substr($_,21,6);
	    $resnum=~s/\s+//g;
	    #print "$resnum $old_resnum $atomtype\n";
	    if($atomtype=~/CA/ && $old_resnum ne $resnum)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	#	print $table{$res};
		$old_resnum=$resnum;
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
   # print "\n";
    return $seq;
}
sub aa321ANY
{
    my $file=shift;
    my $seq="";
    my $old_resnum="whatever";
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 12, 4);
	    my $resnum=substr($_,21,6);
	    $resnum=~s/\s+//g;
	    #print "$resnum $old_resnum $atomtype\n";
	    if($old_resnum ne $resnum)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	#	print $table{$res};
		$old_resnum=$resnum;
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
   # print "\n";
    return $seq;
}

sub aa321_resnum
{
    my $file=shift;
    my $seq="";
    my $old_resnum="whatever";
    my @resnum=();
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM\s/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 12, 4);
	    my $resnum=substr($_,22,5);
#	    $resnum=~s/\s+//g;
	    #print "$resnum $old_resnum $atomtype\n";
	    if($atomtype=~/CA/ && $old_resnum ne $resnum)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	#	print $table{$res};
		push(@resnum,$resnum);
		$old_resnum=$resnum;
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
    if(scalar(@resnum)==0) { #read sequence
	$seq=`grep -v '>' $file`;
	$seq=~s/\n//g;
	my @seq=split(//,$seq);
	for(my $i=0;$i<scalar(@seq);$i++) {
	    push(@resnum,$i+1);
	}
    }

   # print "\n";
    return ($seq,\@resnum);
}

sub parse_TMscore 
{
    my $file=shift;
    my %data=();
    open(TM,$file);
    while(<TM>) {
	if(/^RMSD/)    {
	    my @t=split(/\s+/);
	    $data{'rmsd'}=$t[5];
	}
	if(/^TM-score/)    {
	    my @t=split(/\s+/);
	    $data{'TM'}=$t[2];
	}
	if(/^MaxSub-score/)    {
	    my @t=split(/\s+/);
	    $data{'MX'}=$t[1];
	}
	if(/^GDT-TS/)    {
	    my @t=split(/\s+/);
	    $data{'GDT'}=$t[1];
	}
	if(/^GDT-HA/)    {
	    my @t=split(/\s+/);
	    $data{'GDT_HA'}=$t[1];
	}
	if(/^Length=\s+(\d+)/)    {
	    $data{'GDT_HA'}=$1;
	}
	
    }
    close(TM);
    return(%data);

}
sub aa321_resnumANY
{
    my $file=shift;
    my $seq="";
    my $old_resnum="whatever";
    my @resnum=();
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 12, 4);
	    my $resnum=substr($_,22,5);
#	    $resnum=~s/\s+//g;
	    #print "$resnum $old_resnum $atomtype\n";
	    if($old_resnum ne $resnum)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	#	print $table{$res};
		push(@resnum,$resnum);
		$old_resnum=$resnum;
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
   # print "\n";
    return ($seq,\@resnum);
}

sub aa321CA_legacy
{
    my $file=shift;
    my $seq="";
    open(PDB,"$file");
    while(<PDB>)
    {
	if(/^ATOM/)
	{
	    my $atomno=substr($_, 7, 4);
	    my $atomtype=substr($_, 13, 3);
	    if($atomtype=~/CA/)
	    {
		$res=substr($_,17, 3);
		$seq.=aa321($res);
	    }
	}
	
	last if(/^ENDMDL/);
    }
    close(PDB);
   # print "\n";
    return $seq;
}


sub aa321
{
    my $aa=shift;
    my %aa321=('ALA', 'A',
	       'ARG', 'R',
	       'ASN', 'N',
	       'ASP', 'D',
	       'CYS', 'C',
	       'GLN', 'Q',
	       'GLU', 'E',
	       'GLY', 'G',
	       'HIS', 'H',
	       'ILE', 'I',
	       'LEU', 'L',
	       'LYS', 'K',
	       'MET', 'M',
	       'PHE', 'F',
	       'PRO', 'P',
	       'SER', 'S',
	       'THR', 'T',
	       'TRP', 'W',
	       'TYR', 'Y',
	       'VAL', 'V',
	       'ASX', 'B',
	       'GLX', 'Z',
	       'XXX', 'A',
	       'MSE', 'M',
	       'FME', 'M',
	       'PCA', 'E',
	       '5HP', 'E',
	       'SAC', 'S',
	       'CCS', 'C');
    my $aa1=0;
    $aa1=$aa321{$aa} if(defined($aa321{$aa}));
    return($aa1);

}

sub aa123
{
    my $aa=shift;
    my %aa123=('A','ALA',
	       'R','ARG', 
	       'N','ASN', 
	       'D','ASP',
	       'C','CYS',
	       'Q','GLN', 
	       'E','GLU', 
	       'G','GLY', 
	       'H','HIS',
	       'I','ILE',
	       'L','LEU', 
	       'K','LYS', 
	       'M','MET', 
	       'F','PHE',
	       'P','PRO',
	       'S','SER', 
	       'T','THR', 
	       'W','TRP',
	       'Y','TYR', 
	       'V','VAL',
	       'B','ASX',
	       'Z','GLX');
    my $aa3=0;
    if(defined($aa123{$aa}))
    {
	 $aa3=$aa123{$aa}
    }
    else
    {
	print STDERR "\"$aa\" is not defined!\n"; 
    }
    return($aa3);
}
sub sparse_encode
{
    my $seq=shift;
    #print $seq."\n";
    my @seq=split(//,$seq);
    my $code="";
    my $output="";
    my %aaspe =('A', "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'R', "0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'N', "0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'D', "0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'C', "0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'Q', "0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'E', "0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'G', "0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'H', "0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ",
		'I', "0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ",
		'L', "0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ",
		'K', "0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ",
		'M', "0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ",
		'F', "0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ",
		'P', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ",
		'S', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ",
		'T', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ",
		'W', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ",
		'Y', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ",
		'V', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ",
		'X', "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ",
		);
    
    my $i=0;
    for($i=0;$i<scalar(@seq);$i++)
    {
	$output.=$aaspe{$seq[$i]};
    }
    my @temp=split(/\s+/,$output);
   # print scalar(@temp);
   # print " $i\n";
    if(scalar(@temp)==21*$i)
    {
	chop($output);
	$code=$output;
    }


    return $code;

}



sub ali2segmod
{
    my ($target_ali,$target_id,$template_ali,$template_id)=@_;
    my $return_target="";
    my $return_template="";
    my @target_ali=split(//,$target_ali);
    my @template_ali=split(//,$template_ali);
    my $counter=0;
    while(scalar @target_ali<100)
    {
	push(@target_ali,"-");
	push(@template_ali,"-");
    }
    my $len=scalar @target_ali;#length($target_ali);
    for(my $i=0;$i<scalar @target_ali;$i++)
    {
	$return_target.=" $target_ali[$i]";
	$return_template.=" $template_ali[$i]";
	$counter++;
	#printf("%d\n",length($return_targ% 60);
	$return_target.="\n" if($counter % 30 == 0 && $i!=scalar @target_ali-1); 
	$return_template.="\n" if($counter % 30 == 0 && $i!=scalar @target_ali-1);
    }
    my $str="$len $template_id tit2\n$return_template\n$len $target_id tit2\n$return_target\n";
    #print "$len $template_id tit2\n$return_template\n";
    #print "$len $target_id tit2\n$return_target\n";
    return $str;
}

sub ali2whatif
{
    my ($target_ali,$target_id,$template_ali,$template_id)=@_;
    my $return_target="";
    my $return_template="";
    my @target_ali=split(//,$target_ali);
    my @template_ali=split(//,$template_ali);
    my $len=length($target_ali);
    my $counter=0;
    for(my $i=0;$i<scalar @target_ali;$i++)
    {
	$return_target.="$target_ali[$i]";
	$return_template.="$template_ali[$i]";
	$counter++;
	#printf("%d\n",length($return_targ% 60);
	$return_target.="\n" if($counter % 60 == 0); 
	$return_template.="\n" if($counter % 60 == 0);
    }
    $return_target.="*";
    $return_template.="*";
    my $target_str=">P1;$target_id\ntarget\n$return_target\n";
    my $template_str=">P1;$template_id\ntemplate\n$return_template\n";
    #print "$len $template_id tit2\n$return_template\n";
    #print "$len $target_id tit2\n$return_target\n";
   # print $target_str."\n".$template_str."\n";
    return ($target_str,$template_str);
}

1; #return true

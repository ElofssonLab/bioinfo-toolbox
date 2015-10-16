#!/usr/bin/perl -w

# Note:
# mySCAMPI_run2.pl, without using xslt, which can not handel huge xml files
use File::Basename;
$scampi_dir="//scampi";
$inputFileFormat="fa";

$usage="
Usage: mySCAMPI_run.pl [options] fasta-seq-file
Options:
    --tmpdir <dir>    : set temporary directory
    --not-clean|-nc   : do not clean the temporary dir
    --modhmmopt <str> : add modhmms option
    --outpath  <dir>  : output the result to outpath
    -f|-format        : format of the input file, can be fa, msa, s, prf

Created 2010-08-16, updated 2010-08-23, Nanjiang 
";

sub PrintHelp {
    print $usage;
}

$numArgs = $#ARGV+1;
if($numArgs < 1) {
    &PrintHelp;
    exit;
}

$infile_withpath="";
$outpath="";
$isClean=1;
$tmpdir="";
$modhmms_options="";
if(@ARGV)#{{{
{
    $i = 0;
    while($ARGV[$i])
    {
        if($ARGV[$i] eq "-h" || $ARGV[$i] eq "--help" )
        {
            &PrintHelp;
            exit;
        }
        elsif($ARGV[$i] eq "-tmpdir" || $ARGV[$i] eq "--tmpdir" )
        {   
            $tmpdir = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-outpath" || $ARGV[$i] eq "--outpath" )
        {   
            $outpath = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-f" || $ARGV[$i] eq "--format" )
        {   
            $inputFileFormat = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-modhmmopt" || $ARGV[$i] eq "--modhmmopt" ||
            $ARGV[$i] eq "-extra")
        {   
            $modhmms_options = $modhmms_options . " " . $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "--not-clean" || $ARGV[$i] eq "-nc" )
        {   
            $isClean = 0;
            $i +=1;
        }
        else 
        {
            $infile_withpath = $ARGV[$i];
            $i ++;
        }
    }
}#}}}

$infile=$infile_withpath;
$infile=~ s/.*\///;
$infile_path=dirname($infile_withpath);
$infile_basename=basename($infile_withpath);
if ($outpath eq "")
{
    $outpath=$infile_path;
}
else
{
    system("mkdir -p $outpath");
}

if ($tmpdir eq "")
{
    $tmpdir=`/bin/mktemp -d /tmp/SCAMPI_XXXXXXXXXX`;
    chomp($tmpdir);
}
else
{
    system("/bin/mkdir -p \"$tmpdir\"");
}
print "tmpdir=$tmpdir\n";
# print "infile_path=$infile_path\n";
# print "infile_basename=$infile_basename\n";
#print "isClean=$isClean\n";
#print "modhmms_options=$modhmms_options\n";

mkdir("$tmpdir/sequences");
system("/bin/cp $infile_withpath $tmpdir");
system("/usr/bin/perl $scampi_dir/fasta_total_split.pl $tmpdir/$infile $tmpdir/sequences/seq.");
system("/bin/ls -1 $tmpdir/sequences/|/bin/awk '{print \"$tmpdir/sequences/\" \$1}' > $tmpdir/snf");
# No -g flag (global protein filter); corresponds to what's used in paper
#system("/bin/modhmms_scampi -f fa -s $tmpdir/snf -m $scampi_dir/DGHMM_KR_21.txt -o $tmpdir -r $scampi_dir/replacement_letter_multi.rpl -L --nopostout --nolabels --viterbi -u $modhmms_options > $tmpdir/outfile.xml");
print "/bin/modhmms_scampi -f $inputFileFormat -s $tmpdir/snf -m $scampi_dir/DGHMM_KR_21.txt -o $tmpdir -r $scampi_dir/replacement_letter_multi.rpl -L --viterbi -u $modhmms_options > $tmpdir/outfile.xml\n";
system("/bin/modhmms_scampi -f $inputFileFormat -s $tmpdir/snf -m $scampi_dir/DGHMM_KR_21.txt -o $tmpdir -r $scampi_dir/replacement_letter_multi.rpl -L --viterbi -u  $modhmms_options > $tmpdir/outfile.xml");

system("/misc/casiodata3/bin/my_modhmmxml2res.py $tmpdir/outfile.xml -o $tmpdir/DGHMM_KR_21.hmg.res");
system("/usr/bin/perl $scampi_dir/res2compacttopo.pl $tmpdir/DGHMM_KR_21.hmg.res $tmpdir/sequences > $outpath/$infile_basename.res");
system("/misc/casiodata3/bin/scampiXML2TXT.py -m 2 -i $tmpdir/outfile.xml -aapath $tmpdir/sequences > $outpath/$infile_basename.xml.res ");

if ($isClean == 1)
{
    system("/bin/rm -rf $tmpdir");
}

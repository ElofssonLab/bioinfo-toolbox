#!/usr/bin/perl -w
use File::Basename;
use Cwd qw(realpath);

my $fullpath = realpath($0);
my $programpath=dirname($fullpath);

$scampi_dir="/server/scampi";
$modhmm_dir="/server/modhmm/bin";
$maxNrSeqPerSplit = 500; #because xslt can not handel huge sequences

$usage="
Usage: mySCAMPI_run.pl [options] fasta-seq-file
Options:
    Note: input file is fasta-seq-file (labeled or non-labeled)
    --tmpdir        <dir> : set temporary directory
    --not-clean|-nc       : do not clean the temporary dir
    --modhmmopt     <str> : add modhmms option
    --outpath       <dir> : output the result to outpath
    --maxnrseq      <int> : set the maximum number of sequences in one scampi run, default=$maxNrSeqPerSplit
                          : if the fasta file contains more sequences, it will be splitted

Created 2010-08-16, updated 2015-02-13, Nanjiang 

Examples: 
    ./mySCAMPI_run.pl test/test.fa --outpath out1
";

sub PrintHelp
{
    print $usage;
}
sub SplitSeqFileList#{{{
{
    ($totalSeqFileListFile,$snfFileListFile, $maxNrSeqPerSplit , $tmpdir) = @_;
    open (IN, "<$totalSeqFileListFile") || die "Can not open $totalSeqFileListFile for read\n";
    $cntLine=0;
    $cntSNFFile=1;
    while (<IN>)
    {
        $line=$_;
        chomp($line);
        if ($cntLine == 0)
        {
            open(OUT,">$tmpdir/snf.$cntSNFFile") || die " Can not open $tmpdir/snf.$cntSNFFile for write\n";
        }
        print OUT "$line\n";
        $cntLine++;
        if ($cntLine >= $maxNrSeqPerSplit)
        {
            close(OUT);
            $cntLine=0;
            $cntSNFFile++;
        }
    }
    if ($cntLine != 0)
    {
        close(OUT);   
        $numSNFFile=$cntSNFFile;
    }
    else
    {
        $numSNFFile=$cntSNFFile-1;
    }

    open (OUT, ">$snfFileListFile") || die "Can not open $snfFileListFile for write\n";
    for ($i =1; $i <=$numSNFFile; $i++)
    {
        printf OUT "$tmpdir/snf.%d\n",$i;
    }
    close(OUT);
    
    return $numSNFFile;
}#}}}

$numArgs = $#ARGV+1;
if($numArgs < 1)
{
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

#system("/bin/ls -1 $tmpdir/sequences/|/bin/awk '{print \"$tmpdir/sequences/\" \$1}' > $tmpdir/snf");
system("/bin/ls -1 $tmpdir/sequences/|awk '{print \"$tmpdir/sequences/\" \$1}' > $tmpdir/totalseqfilelist");

$numSNFFile=SplitSeqFileList("$tmpdir/totalseqfilelist", "$tmpdir/snflist", $maxNrSeqPerSplit, $tmpdir);

# No -g flag (global protein filter); corresponds to what's used in paper
#system("/bin/modhmms_scampi -f fa -s $tmpdir/snf -m $scampi_dir/DGHMM_KR_21.txt -o $tmpdir -r $scampi_dir/replacement_letter_multi.rpl -L --nopostout --nolabels --viterbi -u $modhmms_options > $tmpdir/outfile.xml");
system("/bin/cat /dev/null > $outpath/$infile_basename.res");
system("/bin/cat /dev/null > $outpath/$infile_basename.xml.res");
for ($i = 1; $i <= $numSNFFile; $i ++)
{
    $snfFile=sprintf("$tmpdir/snf.%d", $i);
    $tmpXMLFile=sprintf("$tmpdir/outfile.%d.xml", $i);
    $tmpResFile=sprintf("$tmpdir/DGHMM_KR_21.hmg.%d.res", $i);
    print "$modhmm_dir/modhmms_scampi -f fa -s $snfFile -m $scampi_dir/DGHMM_KR_21.txt -o $tmpdir -r $scampi_dir/replacement_letter_multi.rpl   --viterbi -u --labeloddsout --labelllout --labelrevout --nopostout -g --alignlabelout $modhmms_options > $tmpXMLFile\n";
    system("$modhmm_dir/modhmms_scampi -f fa -s $snfFile -m $scampi_dir/DGHMM_KR_21.txt -o $tmpdir -r $scampi_dir/replacement_letter_multi.rpl  --viterbi -u --labeloddsout --labelllout --labelrevout --nopostout -g --alignlabelout  $modhmms_options > $tmpXMLFile");
    system("$modhmm_dir/modhmmxml2res < $tmpXMLFile > $tmpResFile");
    system("/usr/bin/perl $scampi_dir/res2compacttopo.pl $tmpResFile $tmpdir/sequences >> $outpath/$infile_basename.res");
    system("$programpath/scampiXML2TXT.py -m 1 -i $tmpXMLFile -aapath $tmpdir/sequences >> $outpath/$infile_basename.xml.res ");
}
system("$scampi_dir/compacttopo2top.py $outpath/$infile_basename.res > $outpath/$infile_basename.topo" );

if ($isClean == 1)
{
    system("/bin/rm -rf $tmpdir");
}

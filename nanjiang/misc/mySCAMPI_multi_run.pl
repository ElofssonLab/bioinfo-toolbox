#!/usr/bin/perl -w

use strict;
use File::Temp "tempdir";

my $blast_db = "uniprotWG100Filtered";
my $blastdir = $ENV{"BLASTBIN"};
$blastdir="$blastdir/../";
my $outfile ="/dev/stdout";

my $usage="
Usage:  mySCAMPI_multi_run.pl [Options] <fasta-seq-file>
Options:
    Note: input file is fasta-seq-file (labeled or non-labeled)
    --blastdb   <str>  : set the blastdb, default=$blast_db
    --blastdir  <str>  : set the blastdir, default=$blastdir
    --modhmmopt <str>  : add modhmms option
    --outfile   <file> : output the result to outpath, default=$outfile
    --tmpdir    <dir>  : set temporary directory
    --not-clean|-nc    : do not clean the temporary dir

Created 2010-08-27, updated 2010-08-27, Nanjiang 
";

my $numArgs = $#ARGV+1;
if($numArgs < 1)
{
    &PrintHelp;
    exit;
}

# if ( $#ARGV < 3 ) {
#     printf "Usage: $0 <fasta file> <outfile> <BLAST installation directory> <database>\n";
#     exit(1);
# }

my $infile_withpath="";

my $isClean=1;
my $modhmms_options="";

my $bindir = "//scampi-msa";
my $tmpdir = tempdir("/tmp/scampi-msa_tmpdir_XXXXXX");

if(@ARGV)#{{{
{
    my $i = 0;
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
        elsif($ARGV[$i] eq "-o" || $ARGV[$i] eq "--outfile"  || $ARGV[$i] eq "-outfile" )
        {   
            $outfile = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-blastdb" || $ARGV[$i] eq "--blsatdb" )
        {   
            $blast_db = $ARGV[$i+1];
            $i +=2;
        }
        elsif($ARGV[$i] eq "-blastdir" || $ARGV[$i] eq "--blastdir" )
        {   
            $blastdir = $ARGV[$i+1];
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


system("mkdir -p $tmpdir");
my $fastafile=$infile_withpath;

$ENV{'PATH'}="/bin:/usr/bin";
system("$bindir/fix_fasta.pl $fastafile > $tmpdir/query.fa");
#system("/bin/echo $bindir/DGHMM_KR_21_multi.hmg > $bindir/DGHMM_KR_21_multi.txt");

print "BLASTing...\n";
system("$bindir/fa2prfs.sh $tmpdir/query $blastdir $blast_db");
print "BLAST finished!\n";

system("/bin/echo $tmpdir/query.raw.prf > $tmpdir/query.raw.prf.snf");
system("/bin/echo query > $tmpdir/query.pnf");

if(open(IN,"$tmpdir/query.fix")) {
    my $fix_topo = <IN>;
    close(IN);
    $fix_topo =~ /^[iMo\.]+$/ || die;
    relabel_prf_file("$tmpdir/query.prf",$fix_topo);
    relabel_prf_file("$tmpdir/query.raw.prf",$fix_topo);
}

print "tmpdir=$tmpdir\n";
print "Running SCAMPI \n";
system("/bin/modhmms_scampi -f prf -s $tmpdir/query.raw.prf.snf -m $bindir/DGHMM_KR_21_multi.txt -r $bindir/replacement_letter_multi.rpl --nopostout --viterbi -u -L -g $modhmms_options > $tmpdir/scampi_modhmmres.xml");
system("/bin/modhmmxml2res < $tmpdir/scampi_modhmmres.xml > $tmpdir/scampi_modhmmres.res");
system("/data3/bin/swapfile $tmpdir/query.fa $tmpdir/query.raw.prf");
system("/usr/bin/perl $bindir/res2compacttopo.pl $tmpdir/scampi_modhmmres.res $tmpdir > $tmpdir/scampi_modhmmres.compacttopo");
system("$bindir/compacttopo2top.py $tmpdir/scampi_modhmmres.compacttopo > $outfile");
system("/data3/bin/swapfile $tmpdir/query.fa $tmpdir/query.raw.prf");
if ($isClean == 1)
{
    system("/bin/rm -rf $tmpdir");
}


sub relabel_prf_file {#{{{
    my $prffile = shift;
    my $topology = shift;
    if($prffile eq "" || $topology eq "") {
	return;
    }
    else {
	my @topology = split //, $topology;
	
	
	open(PRFFILE, "$prffile")
	    or die "Could not open $prffile";
	
	my $intro = "";
	my @cols = ();
	my $extro = "";
	my $inintro = 1;
	my $inextro = 0;
	my $topopos = 0;
	while(<PRFFILE>) {
	    if($_ =~ "START 1") {
		$intro .= $_;
		$inintro = 0;
	    }
	    if($_ =~ "END 1") {
		$inextro = 1;
	    }
	    if($_ =~ 'ALPHABET' && $inextro <= 0) {
		$intro .= $_;
	    }
	    if($_ =~ 'COL ') {
		chomp;
		$_ =~ s/\s+/ /g;
		my @col = split /\s/, $_;
		$col[$#col - 1] = $topology[$topopos];
		push @cols, [@col];
		$topopos++;
	    }
	    if($inintro > 0) {
		$intro .= $_;
	    }
	    if($inextro > 0) {
		$extro .= $_;
	    }
	    
	}
	close PRFFILE;
	open OUTFILE, ">"."$prffile"
	    or die "Could not open $prffile for writing";
	print OUTFILE "$intro";
	for(my $i = 0; $i <= $#cols; $i++) {
	    my @col = @{$cols[$i]};
	    printf OUTFILE "COL %4d", $i+1;
	    print OUTFILE ":  ";
	    for $b (2 .. $#col ) {
		if($b >= ($#col - 1)) {
		    print OUTFILE "$col[$b]       ";
		}
		else {
		    printf OUTFILE "%5.02f   ", $col[$b];
		}
	    }
	    print OUTFILE "\n";
	}
	print OUTFILE "$extro";
	close OUTFILE;
    }
    return;
}#}}}
sub PrintHelp
{
    print $usage;
}

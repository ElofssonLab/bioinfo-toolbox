use Cwd 'abs_path';
use File::Basename;

my $LWP_loaded=1;
my $HTML_loaded=1;

eval {
    require LWP::UserAgent; 
    LWP::UserAgent->import();
    1;
} or do {
   my $error = $@;
   $LWP_loaded=0

};

eval {
    require HTML::Parser;
    HTML::Parser->import();
    1;
} or do {
   my $error = $@;
   $HTML_loaded=0;
};

my $progname = basename($0);
my $progname_15char = substr $progname, 0, 15;

my $usage ="
usage: $progname -i seqfile [-outpath DIR] [-format html|txt] [-o OUTFILE]

Description: 
    Run topcons on topcons.net

OPTIONS:
    -outpath DIR        Download the result locally to outpath
    -format html|txt    output the result in html or txt format, (default: html)
    -o OUTFILE          output the result to OUTFILE, (default: stdout)

Examples:
    $progname -i test.fa -outpath result -o link.html

Created 2014-10-31, updated 2014-10-31, Nanjiang Shu
";

my $numArgs = $#ARGV+1;
if($numArgs < 1) {
    print "$usage\n";
    exit;
}

my $seqfile = "";
my $outpath = "";
my $format = "html";
my $outfile = "";
if(@ARGV){#{{{
    my $i = 0;
    while($ARGV[$i]) {
        if($ARGV[$i] eq "-h" || $ARGV[$i] eq "--help" ) {
            print "$usage\n";
            exit;
        } elsif($ARGV[$i] eq "-outpath" || $ARGV[$i] eq "--outpath" ) {
            $outpath = File::Spec->rel2abs($ARGV[$i+1]);
            $i += 2;
        } elsif($ARGV[$i] eq "-i" || $ARGV[$i] eq "--i" ) {
            $seqfile = File::Spec->rel2abs($ARGV[$i+1]);
            $i += 2;
        } elsif($ARGV[$i] eq "-o" || $ARGV[$i] eq "--o" ) {
            $outfile = File::Spec->rel2abs($ARGV[$i+1]);
            $i += 2;
        } elsif($ARGV[$i] eq "-format" || $ARGV[$i] eq "--format" ) {
            $format = $ARGV[$i+1];
            $i += 2;
        } else {
            die "wrong argument $ARGV[$i]";
        }
    }
}#}}}

die "$0: $seqfile does not set\n" if ($seqfile eq "");
die "$0: $outpath does not set\n" if ($outpath eq "");
die "$0: unrecognized format \"$format\"\n" if ($format ne "html" || $format ne "txt");


if($LWP_loaded && $HTML_loaded) {
    $topcons="$pdb.topcons";
    $topcons_fa="$pdb.topcons.fa";
    if($overwrite || !-e $topcons) {
        print "topcons...\n";
        my $ua = new LWP::UserAgent;
        my $response = $ua -> post('http://topcons.net',{'sequence' => $seq,'do' => 'Submit',});
        $content=$response->content();
        if($content=~/result\/(.+)\/topcons.txt/) {
            $id=$1;
            print $id."\n";
            my $response = $ua -> post("http://topcons.net/result/$id/topcons.txt");
            open(OUT,">$topcons");
            print OUT $response->content();
            close(OUT);
            #my $parser = new MyParser;
#	my $parsed=$parser->parse($content);
#	print $parsed."\n";
#	$string="sequence=$seq\&do=Submit";
#	$output=`echo "$string" | lynx -post_data http://topcons.net`;
#	    exit;
#	open(OUT,">$topcons");
#	print OUT $output;
#	close(OUT);
        }
#`run_topcons.pl $fasta > $topcons`;
    }

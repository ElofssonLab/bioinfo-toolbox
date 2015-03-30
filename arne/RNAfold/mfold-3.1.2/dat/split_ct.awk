# Split a ct file into individual files.
# Keyword for recognizing a new folding is "="

BEGIN	{
	count = 0;
	if ( ARGC-- != 3 ) {
	  print "usage:";
	  print "nawk -f split.awk file-to-split prefix-for-split-files";
	  exit 1;
	  }
	CurrentFile = ARGV[2];
 	}

$3  ~ /\=/ {
	if ( CurrentFile != ARGV[2] ) close(CurrentFile);
	count++;
	CurrentFile = ARGV[2] "_" count ".ct";
	print "CurrentFile=",CurrentFile;
	}

{ print $0 > CurrentFile; }

END	{ if ( count ) close(CurrentFile); }


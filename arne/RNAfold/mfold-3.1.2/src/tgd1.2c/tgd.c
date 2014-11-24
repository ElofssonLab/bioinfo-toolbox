/* ----------------------------------------------------------
 *
 *	tgd  --a text frontend to the gd graphics library
 *
 *	To do:
 *
 *	o	Improve memory usage, more malloc() less static
 *	
 *	o	Fonts should be dynamically loaded/unloaded
 *
 *	For external documentation see
 *	http://s27w007.pswfs.gov/tgd/
 *
 *	For documentation of the gd library see:
 *	http://siva.cshl.org/gd/gd.html
 *	
 *	Bradley K. Sherman 1995 (bks@netcom.com bks@s27w007.pswfs.gov)
 *
 *	Change History:
 *		Added fonts tiny, mediumbold and giant and
 *		made font names case-insensitive 10/95
 *
 *		Solaris doesn't know about strings.h, deleted
 *		with no adverse effects (I think --bks 10/12/95)
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef NEVER
/* Not available on Solaris, I think it can be discarded */
#include <strings.h>
#endif
#include <ctype.h>
#include <malloc.h>
#include "gd.h"
#include "gdfontl.h"
#include "gdfonts.h"
#include "gdfontg.h"
#include "gdfontmb.h"
#include "gdfontt.h"
#include "tgdsamfont.h"
#include "tgd.h"

static char SCCS[] = "@(#)tgd.c	1.6 10/13/95 tgd";
static char    *Current_file; /* Name of current input file */
static FILE    *Fp_err; /* Where to write errors and debug messages */
static FILE    *Fp_out; /* Where to write standard output */
static TGDFONT  Font[MAXTGDFONT];
static TIM      Tim[MAXTIM]; /* Working images */
static int      Ldebug = 0; /* Debug level set < 0 for no output at all */
static int      N_line = 0;  /* Line in current input file */

COLOUR *colourbyid();
COLOUR *colourbyname();
TIM *timbyname();
static void tgdgo();
void die(char *s );
void tgdputs(char *s );
void tgdputerr(char *s );
void init();
TIM *timopen(char *s );
void tgdstringcenter(int argc, char *argv[] );
void tgdshowvector(int argc, char *argv[] );
void loadfonts();
void unloadfonts();
void timset(TIM *tim );
void tgddebug(int argc,char *argv[] );
void tgdstockimage(int argc,char *argv[] );
void tgdargdestroy(int argc,char *argv[] );

main( argc, argv )
int argc;
char *argv[];
{
	int i;
	FILE *fp;

	if ( argc < 2 ) {
		init();
		Current_file = "STDIN";
		tgdgo( stdin );
	} else {
		for ( i = 1; i < argc; i++ ) {
			if ( NULL == ( fp = fopen( argv[i], "r" ) ) ) {
				perror( argv[i] );
				exit( 1 );
			}
			init();
			Current_file = argv[i];
			tgdgo( fp );
			fclose( fp );
		}
	}
	exit( 0 );
}


/* --
 *	Unrecoverable error, terminate program with
 *	a message.
 */
void die(char *s )
{
	char buf[BUFSIZ];
	sprintf( buf, "tgd error is [%s]. Near line %d in file %s.\n",
		s, N_line, Current_file );
	tgdputerr( buf );
	exit( 1 );
}

/* --
 *	All reporting should go through here.
 *	Setting Ldebug less than 0 will completely
 *	shut program up.
 */
void tgdputs(char *s )

{
	if ( Ldebug >= 0 )
		fprintf( Fp_out, "%s\n", s );
}

/* --
 *	All error reports should go through here.
 */
void tgdputerr(char *s )
{
		fprintf( Fp_err, "%s\n", s );
}

/* --
 *	(Re)Initialize datastructures for each input file.
 */
/*static*/
void init()
{
	int i, j;

	N_line = 0;
	Ldebug = 0;
	Fp_out = stdout;
	Fp_err = stderr;

	unloadfonts(); /* XXX */
	loadfonts(); /* XXX */
	/* Reinit image data structure */
	for (i = 0; i < MAXTIM; i++ )
		if ( Tim[i].im != NULL ) {
			gdImageDestroy( Tim[i].im );
			Tim[i].im = NULL;
			for ( j = 0; j < MAXCOLOUR; j++ )
				strcpy( Tim[i].colour[j].name, "" );
		}
}

/* --
 *  Major loop:
 *	read_command()
 *	tgdcommand()
 */
static void
tgdgo( fp )
FILE *fp;
{
	int argc;
	static char *argv[TGDMAXTOKENS];
	char  buf[TGDMAXLINE];

	/* Get a  line of input and remove newline. */
	while (  NULL != fgets( buf, TGDMAXLINE, fp ) ) {
		++N_line;
		argc = tgdtokenize( buf, argv );
		if ( Ldebug > 90 )
			tgdshowvector( argc, argv );
		tgdcommand( argc, argv );
		tgdargdestroy( argc, argv );
	}
}

/*  ----------------Utility routines follow---------------- */
/* --
 *	If colorname is already allocated for this tim, return its id,
 *	else return -1.
 */
COLOUR *
colourbyname( tim, s )
TIM *tim;
char *s;
{
	int i;
	static COLOUR special[] = {
		{ gdBrushed, "gdBrushed" },
		{ gdBrushed, "gdbrushed" },
		{ gdStyled, "gdStyled" },
		{ gdStyled, "gdstyled" },
		{ gdStyledBrushed, "gdStyledBrushed" },
		{ gdStyledBrushed, "gdstyledbrushed" },
		{ gdTiled, "gdTiled" },
		{ gdTiled, "gdtiled" },
		{ gdTransparent, "gdTransparent" },
		{ gdTransparent, "gdtransparent" },
		{ 91119, "" },  /* changed NULL to "" 10/13/95 --bks */
	};
	for ( i = 0; special[i].id != 91119; i++ )
		if ( strcmp( special[i].name, s ) == 0 )
			return special + i;

	for ( i = 0; i < MAXCOLOUR; i++ )
		if ( strcmp( tim->colour[i].name, s ) == 0 )
			return &(tim->colour[i]);
	return NULL;
}

/* --
 *	Given gd id, return pointer to user color name.
 */
COLOUR *
colourbyid( tim, colourid )
TIM *tim;
int colourid;
{
	int i;
	for ( i = 0; i < MAXCOLOUR; i++ )
		if ( tim->colour[i].id == colourid )
			return &(tim->colour[i]);
	return NULL;
}

/* --
 *	Be sure this TIM is initialized.
 */
TIM *timopen(char *s )
{
	int i, j;
	if ( NULL != timbyname( s ) )
		die( "Duplicate handle in timopen." );
	for ( i = 0; i < MAXTIM; i++ ) {
		if ( Tim[i].name[0] == '\0') {
			strcpy( Tim[i].name, s );
			Tim[i].im = NULL;
			for ( j = 0; j < MAXCOLOUR; j++ ) {
				Tim[i].colour[j].id = -1;
				Tim[i].colour[j].name[0] = '\0';
			}
			return Tim + i;
		}
	}
	die( "too many images" );
}

/* --
 *	Called when we read in an image to synchronize
 *	our conception of the image with the GD library.
 */
void timset(TIM *tim )

{
	int colors;
	int i, r, g, b;
	gdImagePtr im;

	im = tim->im;
	colors = gdImageColorsTotal( im );

	for ( i = 0; i < colors; i++ ) {
		r = gdImageRed( im, i );
		g = gdImageGreen( im, i );
		b = gdImageBlue( im, i );
		sprintf( tim->colour[i].name, "#%02X%02X%02X", r, g, b );
		tim->colour[i].id = i;
	}
}

/* --
 *	Return pointer to TIM specified by name;
 *	create entry for name if it does not yet exist.
 */
TIM *
timbyname( s )
char *s;
{
	int i;
	for ( i = 0; i < MAXTIM; i++ )
		if ( strcmp( Tim[i].name, s ) == 0 )
			return Tim + i;
	return NULL;
}

/* --
 *	Return pointer to TGDFONT with specified name.
 */
TGDFONT *
fontbyname( s )
char *s;
{
	int i;
	char buf[MAXTIMNAME + 1];;
	if ( strlen(s) > MAXTIMNAME )
		die( "Font name too long.");

	strcpy( buf, s );
	for ( i = 0; buf[i] != '\0'; i++ )
		if ( isalpha( buf[i] ) )
			buf[i] = tolower( buf[i] );
	for ( i = 0; i < MAXTGDFONT; i++ )
		if ( strcmp( Font[i].name, buf ) == 0 )
			return Font + i;
	return NULL;
}

/* --
 *	... level
 *	Set verbosity, low value -> quiet
 *	
 *	Setting value to zero (or negative) turns
 *	off all standard output, including that from
 *	printing functions.
 */

void tgddebug(int argc,char *argv[] )
{
	Ldebug = atoi( argv[1] );
}

/* --
 *  center the string on the given x y coords. (see imagestring)
 */
void tgdstringcenter(int argc, char *argv[] )

{
	TIM *tim;
	TGDFONT *font;
	int x, y, len, centerx, centery;
	COLOUR *color;
	char *s;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( font = fontbyname( argv[2] ) ) )
		die( argv[2] );
	x = atoi( argv[3] );
	y = atoi( argv[4] );
	s = argv[5];
	if ( NULL == ( color = colourbyname( tim, argv[6] ) ) )
		die( "Unallocated color for tim." );

	len = strlen(s);
	centerx = x - (font->fo->w * len) / 2;
	centery = y - (font->fo->h) / 2;
	gdImageString( tim->im, font->fo, centerx, centery, s, color->id );
}

/* --
 *	... handle id
 *
 *	Note that "macros" can be built
 *	up out of text commands.
 */
void tgdstockimage(int argc,char *argv[] )

{
	int anotherargc;
	static char *anotherargv[TGDMAXTOKENS]; /* XXX dynamic allocation? */
	char  buf[TGDMAXLINE];
	char **s;

	static char *tgdstock0[] =  {
	"create %s 350 250",
	"colorallocate %s grey    192 192 192",
	"colorallocate %s white   255 255 255",
	"colorallocate %s black     0   0   0",
	"colortransparent %s grey",
	"interlace %s 1",
	NULL,
	};

	static char *tgdstock1[] =  {
	"create %s 350 250",
	"colorallocate %s grey    192 192 192",
	"colorallocate %s white   255 255 255",
	"colorallocate %s black     0   0   0",
	"colorallocate %s red     255   0   0",
	"colorallocate %s orange  255 127   0",
	"colorallocate %s yellow  255 255   0",
	"colorallocate %s green     0 255   0",
	"colorallocate %s blue      0   0 255",
	"colorallocate %s indigo  127   0 255",
	"colorallocate %s violet  255   0 255",
	"colortransparent %s grey",
	"interlace %s 1",
	NULL,
	};

	static char *tgdstock2[] =  {
	"create %s 35 25",
	"colorallocate %s grey    192 192 192",
	"colorallocate %s white   255 255 255",
	"colorallocate %s black     0   0   0",
	"colorallocate %s red     255   0   0",
	"colorallocate %s orange  255 127   0",
	"colorallocate %s yellow  255 255   0",
	"colorallocate %s green     0 255   0",
	"colorallocate %s blue      0   0 255",
	"colorallocate %s indigo  127   0 255",
	"colorallocate %s violet  255   0 255",
	"colortransparent %s grey",
	"interlace %s 1",
	NULL,
	};

	switch( atoi( argv[2] ) ) {
	case 0 : s = tgdstock0; break;
	case 1 : s = tgdstock1; break;
	case 2 : s = tgdstock2; break;
	default : die( "no such stock image id" );
	}

	for ( ; *s != NULL; s++ ) {
		sprintf( buf, *s, argv[1] );
		anotherargc = tgdtokenize( buf, anotherargv );
		tgdcommand( anotherargc, anotherargv ); /* Execute command */
		tgdargdestroy( anotherargc, anotherargv ); /* XXX */
	}
}

/* --
 *	Release argv strings 
 */
void tgdargdestroy(int argc,char *argv[] )

{
	int i;
	for( i = 0; i < argc; i++ ) {
		if (argv[i] != NULL )
			free( argv[i] ); /* XXX */
		argv[i] = NULL;
	}

}

/* --
 *	Report on the current comand.
 */
void tgdshowvector(int argc, char *argv[] )

{
	int i;
	char buf[TGDMAXLINE];

	*buf = '\0';
	for ( i = 0; i < argc; i++ )
		sprintf( buf + strlen( buf ), "<%s>", argv[i] );
	tgdputerr( buf );
}

/* --
 *	Get and tokenize one line of input
 *		tokens in input line are whitespace separated
 *		Double quotes surround single tokens.
 *
 *	Returns  number of tokens found (zero should be okay)
 */
int
tgdtokenize( line, argv )
char *line;
char *argv[];
{
	int   count;
	char *end, *s, *t, *u;
	char  tok[TGDMAXLINE];

	/* Trailing newlines? */
	end = line + strlen( line ) - 1;
	if ( *end == '\n' )
		*end-- = '\0';

	/* Trailing whitespace? */
	while ( isspace( *end ) )
		*end-- = '\0';

	/* Break input into tokens */
	count = 0;
	s = line;
	do {
		while( isspace( *s ) ) ++s; /* eliminate whitespace */
		u = tok;
		if ( *s == DOUBLE_QUOTE ) {
			++s;
			while( *s  &&  *s != DOUBLE_QUOTE )
				*u++ = *s++;
			if ( *s == DOUBLE_QUOTE )
				++s;
			else
				die( "Unterminated string" );
			/* XXX Should allow escape with backslash */
		} else {
			while ( *s  &&  ! isspace( *s ) )
				*u++ = *s++;
		}
		*u = '\0'; /* tok finally contains token, terminate */
		t = malloc( strlen( tok ) + 1 ); /* malloc space for token */
		strcpy( t, tok ); /* copy token into dynamic space */
		argv[count++] = t; /* Load token into argument vector */
	} while( s <= end );
	return count;  /* Return the count of tokens found. */
}

/* _________________________________________________________________ */
/* -----------Below is font stuff; should be separate module-------- */

/* --
 *	Eventually should be dynamic loader
 *
 *	Fonts could be compressed, etc.
 *
 *	Font names are case-insenstive to user, so make names
 *	all lower case here (cf. fontbyname() ).
 */
void loadfonts()
{
	/* Font kludge XXX */
	strcpy( Font[0].name, "gdfontsmall" );
	Font[0].fo = gdFontSmall;
	strcpy( Font[1].name, "gdfontlarge" );
	Font[1].fo = gdFontLarge;
	strcpy( Font[2].name, "gdfont10x20" );
	Font[2].fo = gdFont10x20;
	strcpy( Font[3].name, "gdfont12x24" );
	Font[3].fo = gdFont12x24;
	strcpy( Font[4].name, "gdfont5x8" );
	Font[4].fo = gdFont5x8;
	strcpy( Font[5].name, "gdfont6x10" );
	Font[5].fo = gdFont6x10;
	strcpy( Font[6].name, "gdfont6x12" );
	Font[6].fo = gdFont6x12;
	strcpy( Font[7].name, "gdfont6x13" );
	Font[7].fo = gdFont6x13;
	strcpy( Font[8].name, "gdfont6x13bold" );
	Font[8].fo = gdFont6x13bold;
	strcpy( Font[9].name, "gdfont6x9" );
	Font[9].fo = gdFont6x9;
	strcpy( Font[10].name, "gdfont7x13" );
	Font[10].fo = gdFont7x13;
	strcpy( Font[11].name, "gdfont7x13bold" );
	Font[11].fo = gdFont7x13bold;
	strcpy( Font[12].name, "gdfont7x14" );
	Font[12].fo = gdFont7x14;
	strcpy( Font[13].name, "gdfont8x13" );
	Font[13].fo = gdFont8x13;
	strcpy( Font[14].name, "gdfont8x13bold" );
	Font[14].fo = gdFont8x13bold;
	strcpy( Font[15].name, "gdfont8x16" );
	Font[15].fo = gdFont8x16;
	strcpy( Font[16].name, "gdfont9x15" );
	Font[16].fo = gdFont9x15;
	strcpy( Font[17].name, "gdfont9x15bold" );
	Font[17].fo = gdFont9x15bold;
	/* Added 10/95 --bks */
	strcpy( Font[18].name, "gdfontgiant" );
	Font[18].fo = gdFontGiant;
	strcpy( Font[19].name, "gdfontmediumbold" );
	Font[19].fo = gdFontMediumBold;
	strcpy( Font[20].name, "gdfonttiny" );
	Font[20].fo = gdFontTiny;

}

/* --
 *	Eventually will dump fonts out of memory
 */
void unloadfonts()
{
	/* XXX */
}

/* ---------------------------------------------------------- */

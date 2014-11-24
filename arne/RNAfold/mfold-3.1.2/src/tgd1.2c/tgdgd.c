/* ----------------------------------------------------------
 *
 *	tgdgd.c
 *
 *	This file contains the routines that interface
 *	one-to-one with the gd library routines. 
 *
 *	tgd  --a text frontend to the gd library
 *
 *	For documentation of the gd library see:
 *	http://www.boutell.com/gd/
 *	
 *	Bradley K. Sherman 1995 (bks@netcom.com bks@s27w007.pswfs.gov)
 *
 *	Thanks to Ronan Waide for polygon/filledpolygon fixes 1996
 */

#include <stdio.h>
#include "gd.h"
#include "tgd.h"

static char SCCS[] = "@(#)tgdgd.c	1.4 10/7/96 tgd";

COLOUR *colourbyid();
COLOUR *colourbyname();
TIM *timbyname();
TIM *timopen();
TGDFONT *fontbyname();

/* ----------------------------------------------------------------- */

/* --
 *	... handle  cx cy w h start_degrees end_degrees color
 */
imagearc( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int cx, cy, w, h, s, e;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	cx = atoi( argv[2] );
	cy = atoi( argv[3] );
	w = atoi( argv[4] );
	h = atoi( argv[5] );
	s = atoi( argv[6] );
	e = atoi( argv[7] );
	if ( NULL == ( color = colourbyname( tim, argv[8] ) ) )
		die( argv[6] );
	gdImageArc( tim->im, cx, cy, w, h, s, e, color->id );
	return;
}

/* --
 *  handle color
 */
imageblue( argc, argv )
int argc;
char *argv[];
{

	TIM *tim;
	COLOUR *color;
	char buf[100];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( color = colourbyname( tim, argv[2] ) ) )
		die( argv[2] );
	sprintf( buf, "%d", gdImageBlue( tim->im, color->id ) );
	tgdputs( buf );
}

/* --
 *	... handle x y
 */
void
imageboundssafe( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	char x, y;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x = atoi( argv[2] );
	y = atoi( argv[3] );
	if ( gdImageBoundsSafe( tim->im, x, y ) )
		tgdputs("1");
	else
		tgdputs("0");
}

/* --
 *	... handle font x y c color
 *     (not sure whether char should be an integer, string for now.)
 */
imagechar( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	TGDFONT *font;
	int x, y, c;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( font = fontbyname( argv[2] ) ) )
		die( argv[2] );
	x = atoi( argv[3] );
	y = atoi( argv[4] );
	c = argv[5][0]; /* 1-96?*/
	if ( NULL == ( color = colourbyname( tim, argv[6] ) ) )
		die( argv[6] );
	gdImageChar( tim->im, font->fo, x, y, c, color->id );
}

/* --
 *	... handle font x y c color
 */
imagecharup( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	TGDFONT *font;
	int x, y, c;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( font = fontbyname( argv[2] ) ) )
		die( argv[2] );
	x = atoi( argv[3] );
	y = atoi( argv[4] );
	c = argv[5][0]; /* 1-96?*/
	if ( NULL == ( color = colourbyname( tim, argv[6] ) ) )
		die( argv[6] );
	gdImageCharUp( tim->im, font->fo, x, y, c, color->id );
}

/* --
 *  ... handle colorname red green blue
 */
imagecolorallocate( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int r, g, b;
	int i, colorid;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	r = atoi( argv[3] );
	g = atoi( argv[4] );
	b = atoi( argv[5] );
        colorid = gdImageColorAllocate( tim->im, r, g, b );
	if ( colorid == -1 )
		die( "Couldn't allocate color." );
	for ( i = 0; i < MAXCOLOUR; i++ )
		if ( tim->colour[i].id == -1 ) {
			tim->colour[i].id = colorid;
			strcpy( tim->colour[i].name, argv[2] );
			return;
		}
	die( "can't find space for color" );
}

/* --
 *	Prints closest user colorname to specified rgb.
 *	... handle r g b
 */
imagecolorclosest( argc, argv )
int argc;
char *argv[];
{
	int r, g, b;
	COLOUR *color;
	int  close;
	TIM *tim;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	r = atoi( argv[2] );
	g = atoi( argv[3] );
	b = atoi( argv[4] );

	close = gdImageColorClosest( tim->im, r, g, b );
	if ( NULL == ( color = colourbyid( tim, close ) ) )
		die( "imagecolorclosest" );
	tgdputs( color->name );
}

/* --
 *	... handle color
 */
void
imagecolordeallocate( argc, argv )
int argc;
char *argv[];
{
	COLOUR *color;
	TIM *tim;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( color = colourbyname( tim, argv[2] ) ) )
		die( argv[2] );
	gdImageColorDeallocate( tim->im, color->id );
	color->id = -1;
	color->name[0] = '\0';
}

/* --
 *	... handle r g b
 *	Prints color name matching r g b,
 *	or "-1" if no match.
 */
imagecolorexact( argc, argv )
int argc;
char *argv[];
{
	int r, g, b;
	COLOUR *color;
	int  exact;
	TIM *tim;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	r = atoi( argv[2] );
	g = atoi( argv[3] );
	b = atoi( argv[4] );

	exact = gdImageColorExact( tim->im, r, g, b );
	if ( exact == -1 )
	{
		tgdputs( "-1" );
	} else {
		if ( NULL == ( color = colourbyid( tim, exact ) ) )
			die( "imagecolorexact " );
		tgdputs( color->name );
	}
}

/* --
 *	... handle color
 */
imagecolortransparent( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( color = colourbyname( tim, argv[2] ) ) )
		die( argv[2] );
	gdImageColorTransparent( tim->im, color->id );
}

/* --
 *	... dest_handle source_handle dest_x dext_y src_x src_y w h
 */
void
imagecopy( argc, argv )
int argc;
char *argv[];
{
	TIM *src, *dst;
	int dx, dy, sx, sy, w, h ;
	if ( NULL == ( dst = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( src = timbyname( argv[2] ) ) )
		die( argv[2] );
	dx = atoi( argv[3] );
	dy = atoi( argv[4] );
	sx = atoi( argv[5] );
	sy = atoi( argv[6] );
	w  = atoi( argv[7] );
	h  = atoi( argv[8] );
	gdImageCopy( dst->im, src->im, dx, dy, sx, sy, w, h );
}

/* --
 *	... dst src dest_x dext_y src_x src_y dst_w dst_h src w src h
 */
void
imagecopyresized( argc, argv )
int argc;
char *argv[];
{
	TIM *src, *dst;
	int dx, dy, sx, sy, dw, dh, sw, sh;
	if ( NULL == ( dst = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( src = timbyname( argv[2] ) ) )
		die( argv[2] );
	dx = atoi( argv[3] );
	dy = atoi( argv[4] );
	sx = atoi( argv[5] );
	sy = atoi( argv[6] );
	dw = atoi( argv[7] );
	dh = atoi( argv[8] );
	sw = atoi( argv[9] );
	sh = atoi( argv[10] );
	gdImageCopyResized( dst->im, src->im, dx, dy, sx, sy, dw, dh, sw, sh );
}

/* --
 *	... handle x y
 */
imagecreate( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	tim = timopen( argv[1] );
	tim->im = gdImageCreate( atoi(argv[2]), atoi(argv[3]) );
	if ( NULL == tim->im  ) die( "gdImageCreate" );
}

/* --
 *   ... handle filename
 */
imagecreatefromgd( argc, argv )
int argc;
char *argv[];
{

	FILE *fp;
	TIM  *tim;
	tim = timopen( argv[1] );
	if ( NULL == ( fp = fopen( argv[2], "rb" ) ) ) {
		perror( argv[1] );
		die( "Can't open gd file." );
	}
	if ( NULL == ( tim->im = gdImageCreateFromGd( fp ) ) )
		die( "gdImageCreateFromGd" );
	fclose( fp );
	timset( tim );
	return;
}

/* --
 *   ... handle filename
 */
imagecreatefromgif( argc, argv )
int argc;
char *argv[];
{
	FILE *fp;
	TIM  *tim;
	tim = timopen( argv[1] );
	if ( NULL == ( fp = fopen( argv[2], "rb" ) ) ) {
		perror( argv[1] );
		die( "Can't open gif file." );
	}
	if ( NULL == ( tim->im = gdImageCreateFromGif( fp ) ) )
		die( "gdImageCreateFromGif" );
	fclose( fp );
	timset( tim );
	return;
}

/* --
 *	... handle filename.xbm
 *	XXX Not working?
 */
imagecreatefromxbm( argc, argv )
int argc;
char *argv[];
{
	FILE *fp;
	TIM  *tim;
	tim = timopen( argv[1] );
	if ( NULL == ( fp = fopen( argv[2], "rb" ) ) ) {
		perror( argv[1] );
		die( "Can't open xbm file." );
	}
	if ( NULL == ( tim->im = gdImageCreateFromXbm( fp ) ) )
		die( "gdImageCreateFromXbm" );
	fclose( fp );
	timset( tim );
	return;
}

/* --
 *	... handle
 */
imagecolorstotal( argc, argv )
int argc;
char *argv[];
{
	TIM  *tim;
	char buf[BUFSIZ];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	sprintf( buf, "%d", gdImageColorsTotal( tim->im ) );
	tgdputs( buf );
}

/* --
 * supposedly defunct in 1.1.1 kept for backwards compatibility
 */
imagedashedline( argc, argv )
int argc;
char *argv[];
{

	TIM *tim;
	int x1, y1, x2, y2;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x1 = atoi( argv[2] );
	y1 = atoi( argv[3] );
	x2 = atoi( argv[4] );
	y2 = atoi( argv[5] );
	color = colourbyname( tim, argv[6] );
	if ( color == NULL )
		die( argv[6] );
	gdImageDashedLine( tim->im, x1, y1, x2, y2, color->id );
}

/* --
 *	... handle
 */
imagedestroy( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int i;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( "No such tim." );
	gdImageDestroy( tim->im );
	tim->im = NULL;
	tim->name[0] = '\0';
	for ( i = 0; i < MAXCOLOUR; i++ ) {
		tim->colour[i].id = -1;
		tim->colour[i].name[0] = '\0';
	}
}

/* --
 *	... handle x, y, color
 */
void
imagefill( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int x, y;
	COLOUR *color;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x = atoi( argv[2] );
	y = atoi( argv[3] );
	if ( NULL == ( color = colourbyname( tim, argv[4] ) ) )
		die( argv[4] );
	gdImageFill( tim->im, x, y, color->id );
}

/* --
 * ... handle x, y, border, color
 */
void
imagefilltoborder( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int x, y;
	COLOUR *border, *color;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x = atoi( argv[2] );
	y = atoi( argv[3] );
	if ( NULL == ( border = colourbyname( tim, argv[4] ) ) )
		die( argv[4] );
	if ( NULL == ( color = colourbyname( tim, argv[5] ) ) )
		die( argv[5] );
	gdImageFillToBorder( tim->im, x, y, border->id, color->id );
}

/* --
 *	... handle x1, y1, x2, y2, color
 */
imagefilledrectangle( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int x1, y1, x2, y2;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x1 = atoi( argv[2] );
	y1 = atoi( argv[3] );
	x2 = atoi( argv[4] );
	y2 = atoi( argv[5] );
	if ( NULL == ( color = colourbyname( tim, argv[6] ) ) )
		die( argv[6] );
	gdImageFilledRectangle( tim->im, x1, y1, x2, y2, color->id );
}

/* --
 * 	... handle filename
 */
imagegd( argc, argv )
int argc;
char *argv[];
{
	FILE *fp;
	TIM *tim;
	int close_fp = 1;
	if ( 0 == (strcmp( argv[2], "stdout") ) ) {
		fp = stdout;
		close_fp = 0;
	}
	else if ( NULL == ( fp = fopen( argv[2], "wb" ) ) ) {
		perror( argv[2] );
		die( "Can't open output file." );
	}
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	gdImageGd( tim->im, fp );
	if ( close_fp )
		fclose( fp );
}

/* --
 *	... handle
 *	print 1 if interlaced 0 otherwise
 */
imagegetinterlaced( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( gdImageGetInterlaced( tim->im ) )
		tgdputs( "1" );
	else
		tgdputs( "0" );
}

/* --
 * ... handle x y
 */
imagegetpixel( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int x, y;
	int cdex;
	COLOUR *color;
	char buf[100];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x = atoi( argv[2] );
	y = atoi( argv[3] );
	/* Check boundaries here? */
	cdex = gdImageGetPixel( tim->im, x, y );
	color = colourbyid( tim, cdex );
	if ( color == NULL )
		die( "bad pixel, out of bounds?" );
	tgdputs( color->name );
}

/* --
 * Print user name for color which is transparent for this image.
 */
imagegettransparent( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	COLOUR *color;
	int gdcolorid;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	gdcolorid = gdImageGetTransparent( tim->im );
	if ( gdcolorid == -1 ) {
		tgdputs( "-1" );
	}else {
		color = colourbyid( tim, gdcolorid );
		if ( color == NULL )
			die( "color mixup" );
		tgdputs( color->name );
	}
}

/*
 *  imagegif handle filename
 */
imagegif( argc, argv )
int argc;
char *argv[];
{
	FILE *fp;
	TIM *tim;
	int close_fp = 1;
	if ( 0 == (strcmp( argv[2], "stdout" ) ) ) {
		fp = stdout;
		close_fp = 0;
	}
	else if ( NULL == ( fp = fopen( argv[2], "wb" ) ) ) {
		perror( argv[2] );
		die( "Can't open output file." );
	}
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	gdImageGif( tim->im, fp );
	if ( close_fp )
		fclose( fp );
}

/* --
 * 	...  handle color
 */
imagegreen( argc, argv )
int argc;
char *argv[];
{

	TIM *tim;
	COLOUR *color;
	char buf[100];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( color = colourbyname( tim, argv[2] ) ) )
		die( argv[2] );
	sprintf( buf, "%d", gdImageGreen( tim->im, color->id ) );
	tgdputs( buf );
}

/* --
 *	... handle on/off
 *	1 = interlace, 0 = don't
 *  
 */
void
imageinterlace( argc, argv )
int argc;
char *argv[];
{

	TIM *tim;
	int on_off = 0;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( atoi( argv[2] ) != 0 )
		on_off = 1;
	gdImageInterlace( tim->im, on_off );
}

/* --
 *	... handle x1 y1 x2 y2 color
 */
imageline( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int x1, y1, x2, y2;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x1 = atoi( argv[2] );
	y1 = atoi( argv[3] );
	x2 = atoi( argv[4] );
	y2 = atoi( argv[5] );
	color = colourbyname( tim, argv[6] );
	if ( color == NULL )
		die( argv[6] );
	gdImageLine( tim->im, x1, y1, x2, y2, color->id );
}

/* --
 *	handle x1 y1 ... xn yn n color
 */
imagefilledpolygon( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int npoints, i;
	COLOUR *color;
	gdPoint points[256];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	npoints = atoi( argv[argc - 2] );
	if ( npoints > 256  ||  (npoints + npoints != argc - 4)  )
		die( "Too many points or wrongly specifed");
	for ( i = 0; i < npoints; i++ ) {
		points[i].x = atoi( argv[i + i + 2] );
		points[i].y = atoi( argv[i + i + 3] );
	}
	if ( NULL == ( color = colourbyname( tim, argv[argc - 1] ) ) )
		die( argv[argc - 1] );
	gdImageFilledPolygon( tim->im, points, npoints, color->id );
}

/* --
 *  ... handle point-array, npoints, color
 */
void
imagepolygon( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int npoints, i;
	COLOUR *color;
	gdPoint points[256];

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	npoints = atoi( argv[argc - 2] );
	if ( npoints + npoints != argc - 4  )
		die( "Mistake in points specification");
	for ( i = 0; i < npoints; i++ ) {
		points[i].x = atoi( argv[i + i + 2] );
		points[i].y = atoi( argv[i + i + 3] );
	}
	if ( NULL == ( color = colourbyname( tim, argv[argc - 1] ) ) )
		die( argv[argc - 1] );
	gdImagePolygon( tim->im, points, npoints, color->id );
}

/* --
 *	... handle x1 y1 x2 y2 color
 */
void
imagerectangle( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int x1, y1, x2, y2;
	COLOUR *color;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x1 = atoi( argv[2] );
	y1 = atoi( argv[3] );
	x2 = atoi( argv[4] );
	y2 = atoi( argv[5] );
	if ( NULL == ( color = colourbyname( tim, argv[6] ) ) )
		die( argv[6] );
	gdImageRectangle( tim->im, x1, y1, x2, y2, color->id );
}

/* --
 *	... handle color
 */
imagered( argc, argv )
int argc;
char *argv[];
{

	TIM *tim;
	COLOUR *color;
	char buf[100];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( color = colourbyname( tim, argv[2] ) ) )
		die( argv[2] );
	sprintf( buf, "%d", gdImageRed( tim->im, color->id ) );
	tgdputs( buf );
}

/* --
 *	... handle
 */
imagesx( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int width;
	char buf[100];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	width = gdImageSX( tim->im );
	sprintf( buf, "%d", width );
	tgdputs( buf );
}

/* --
 *	... handle
 */
imagesy( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	int height;
	char buf[100];
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	height = gdImageSY( tim->im );
	sprintf( buf, "%d", height );
	tgdputs( buf );
}

/* --
 * XXX
 * ... handle brush_handle
 * Set a special color as the brush
 */
imagesetbrush( argc, argv )
int argc;
char *argv[];
{
	TIM *tim, *brush;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( brush = timbyname( argv[2] ) ) )
		die( argv[2] );
	gdImageSetBrush( tim->im, brush->im );
}

/* --
 *  ... handle x y color
 */
void
imagesetpixel( argc, argv )
int argc;
char *argv[];
{
	int x, y;
	COLOUR *color;
	TIM *tim;

	if ( argc < 5 )
		die( "imagesetpixel args" );
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	x = atoi( argv[2] );
	y = atoi( argv[3] );
	if ( NULL == ( color = colourbyname( tim, argv[4] ) ) )
		die( "Unallocated color for tim." );

	gdImageSetPixel( tim->im, x, y, color->id );
}

/* --
 *	... handle style stylelength
 *	XXX
 */
void
imagesetstyle( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	COLOUR *color;
	int style[256];
	int n, i;

	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	n = atoi( argv[argc - 1] );
	if ( argc != n + 3 )
		die( "parameter disagreement" );
	if ( n > 256 )
		die( "too many colors" );
	if ( strcmp( argv[2], "1" )  ||  strcmp( argv[2], "0" ) ) {
		for ( i = 0; i < n; i++ ) {
			if ( strcmp( argv[i+2], "0" ) == 0 )
				style[i] = 0;
			else if ( strcmp( argv[i+2], "1" ) == 0 )
				style[i] = 1;
			else
				die( argv[i+2] );
		}
	} else {
		for ( i = 0; i < n; i++ ) {
			if ( NULL == ( color = colourbyname( tim, argv[i+2])))
				die( argv[i+2] );
			style[i] = color->id;
		}
	}
	gdImageSetStyle( tim->im, style, n );
}

/* --
 * 	handle tile_handle
 */
imagesettile( argc, argv )
int argc;
char *argv[];
{

	TIM *tim, *tile;
	if ( NULL == ( tim = timbyname( argv[1] ) ) )
		die( argv[1] );
	if ( NULL == ( tile = timbyname( argv[2] ) ) )
		die( argv[2] );
	gdImageSetTile( tim->im, tile->im );
}

/* --
 *	... handle font x y s color
 */
imagestring( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	TGDFONT *font;
	int x, y;
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
	gdImageString( tim->im, font->fo, x, y, s, color->id );
}

/* --
 *	... handle font x y s color
 */
imagestringup( argc, argv )
int argc;
char *argv[];
{
	TIM *tim;
	TGDFONT *font;
	int x, y;
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
	gdImageStringUp( tim->im, font->fo, x, y, s, color->id );
}
/* ---------end of file tgdgd.c------------------------------------- */

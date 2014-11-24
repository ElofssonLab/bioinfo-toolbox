/* ----------------------------------------------------------
 *
 *	tgdcommand.c
 *
 *	tgd  --a text frontend to the gd library
 *	
 *	Bradley K. Sherman 1995 (bks@netcom.com bks@s27w007.pswfs.gov)
 */
 
#include <stdio.h>
 
static char SCCS[] = "@(#)tgdcommand.c	1.4 10/5/95 tgd";

#define BRUSHED  0
#define DASHSIZE  1
#define FONT  2
#define FONTPTR  3
#define IMAGE  4
#define IMAGEARC  5
#define IMAGEBLUE  6
#define IMAGEBOUNDSSAFE  7
#define IMAGECHAR  8
#define IMAGECHARUP  9
#define IMAGECOLORALLOCATE  10
#define IMAGECOLORCLOSEST  11
#define IMAGECOLORDEALLOCATE  12
#define IMAGECOLOREXACT  13
#define IMAGECOLORTRANSPARENT  14
#define IMAGECOPY  15
#define IMAGECOPYRESIZED 16
#define IMAGECREATE  17
#define IMAGECREATEFROMGD  18
#define IMAGECREATEFROMGIF  19
#define IMAGECREATEFROMXBM  20
#define IMAGEDASHEDLINE  21
#define IMAGEDESTROY  22
#define IMAGEFILL  23
#define IMAGEFILLTOBORDER  24
#define IMAGEFILLEDRECTANGLE  25
#define IMAGEGD  26
#define IMAGEGETINTERLACED 27
#define IMAGEGETPIXEL  28
#define IMAGEGETTRANSPARENT  29
#define IMAGEGIF  30
#define IMAGEGREEN  31
#define IMAGEINTERLACE  32
#define IMAGELINE  33
#define IMAGEFILLEDPOLYGON  34
#define IMAGEPOLYGON  35
#define IMAGEPTR  36
#define IMAGERECTANGLE  37
#define IMAGERED  38
#define IMAGESETBRUSH  39
#define IMAGESETPIXEL 40
#define IMAGESETSTYLE  41
#define IMAGESETTILE  42
#define IMAGESTRING  43
#define IMAGESTRINGUP  44
#define MAXCOLORS  45
#define POINT  46
#define STYLED  47
#define STYLEDBRUSHED  48
#define TILED  49
#define TRANSPARENT 50
#define TGDCOMMENT 51
#define TGDDEBUG 52
#define TGDBLANKLINE 53
#define TGDSTOCKIMAGE 54
#define TGDCHECKSYNTAX 55  /* Handled differently from other commands */
#define IMAGESX 56
#define IMAGESY 57
#define IMAGECOLORSTOTAL 58
#define TGDSTRINGCENTER 59

typedef struct cmd {
	char *name;  /* Name supplied by user */
	int   id;  /* Internal identifier */
	char *argcheck;	/* To see if supplied arguments are correct */
} CMD; /* ASCII name id of command */

CMD Command[] = {
    { "arc", IMAGEARC, "ssdddddds" },
    { "blue", IMAGEBLUE, "sss" },
    { "boundssafe", IMAGEBOUNDSSAFE, "ssdd" },
    { "char", IMAGECHAR, "sssddss" }, /* XXX d or s */
    { "charup", IMAGECHARUP, "sssddss" },
    { "colorallocate", IMAGECOLORALLOCATE, "sssddd" },
    { "colorclosest", IMAGECOLORCLOSEST, "ssddd" },
    { "colordeallocate", IMAGECOLORDEALLOCATE, "sss" },
    { "colorexact", IMAGECOLOREXACT, "ssddd" },
    { "colortransparent", IMAGECOLORTRANSPARENT, "sss" },
    { "copy", IMAGECOPY, "sssdddddd" },
    { "copyresized", IMAGECOPYRESIZED, "sssdddddddd" },
    { "create", IMAGECREATE, "ssdd" },
    { "createfromgd", IMAGECREATEFROMGD, "sss" },
    { "createfromgif", IMAGECREATEFROMGIF, "sss" },
    { "createfromxbm", IMAGECREATEFROMXBM, "sss" },
    { "colorstotal", IMAGECOLORSTOTAL, "ss" },
    { "dashedline", IMAGEDASHEDLINE, "ssdddds" },
    { "destroy", IMAGEDESTROY, "ss" },
    { "fill", IMAGEFILL, "ssdds" },
    { "filltoborder", IMAGEFILLTOBORDER, "ssddss" },
    { "filledrectangle", IMAGEFILLEDRECTANGLE, "ssdddds" },
    { "gd", IMAGEGD, "sss" },
    { "getinterlaced", IMAGEGETINTERLACED, "ss" },
    { "getpixel", IMAGEGETPIXEL, "ssdd" },
    { "gettransparent", IMAGEGETTRANSPARENT, "ss" },
    { "gif", IMAGEGIF, "sss" },
    { "green", IMAGEGREEN, "sss" },
    { "interlace", IMAGEINTERLACE, "ssd" },
    { "line", IMAGELINE, "ssdddds" },
    { "filledpolygon", IMAGEFILLEDPOLYGON, "ssddddddd*" },
    { "polygon", IMAGEPOLYGON, "ssddddddd*" },
    { "rectangle", IMAGERECTANGLE, "ssdddds" },
    { "red", IMAGERED, "sss" },
    { "sx", IMAGESX, "ss" },
    { "sy", IMAGESY, "ss" },
    { "setbrush", IMAGESETBRUSH, "sss" },
    { "setpixel", IMAGESETPIXEL, "ssdds" },
    { "setstyle", IMAGESETSTYLE, "ss*" }, /* XXX working? */
    { "settile", IMAGESETTILE, "sss" },
    { "string", IMAGESTRING, "sssddss" },
    { "stringup", IMAGESTRINGUP, "sssddss" },
    { "debug", TGDDEBUG, "sd" },
    { "stockimage", TGDSTOCKIMAGE, "ssd" },
    { "stringcenter", TGDSTRINGCENTER, "sssddss" },
    { "checksyntax", TGDCHECKSYNTAX, "s*" },
    { NULL, -1, NULL }  /* Sentinel. Must be last entry in this array */
};

static CMD *cmdbyid();

/* -------------------------------------------------------------------------- */

/* --
 *	First token on each input line should be a command
 */
static int
which_command( argc, argv )
int argc;
char *argv[];
{
	char buf[BUFSIZ];
	CMD *cmd;
	int i;
	char *s;
	char *t;

	if ( argc == 0 )
		return TGDBLANKLINE;

	s = argv[0];

	if ( (i = strlen(s)) >= BUFSIZ )
		die( "command name too long" );

	if ( i == 0 )
		return TGDBLANKLINE;

	if ( *s == '#' )
		return TGDCOMMENT;

	for( t = buf; *s; s++, t++ ) /* Map to lowercase */
		*t = tolower( *s );
	*t = '\0';

	for ( cmd = Command; cmd->name != NULL; cmd++ )
		if ( strcmp( buf, cmd->name ) == 0 )
			return cmd->id;

	strcat( buf, " is not a known command" );
	die( buf );
}

/* --
 *	Return CMD struct matching id
 */
static CMD *
cmdbyid( id )
int id;
{
	CMD *cmd;
	for ( cmd = Command; cmd->name != NULL; cmd++ )
		if ( cmd->id == id )
			return cmd;
	die( "command not found" );  /* XXX return NULL ? */
}

/* --
 *	Check arguments for type, number
 *	d = integer, s = string, * = any number of following args
 *	zeroth character is always s for the command itself 
 *
 *	Returns 0 on Success, -1 on failure.
 */
static int
tgdargcheck( cmd_id, argc, argv )
int cmd_id;
int argc;
char *argv[];
{
	char *s, *t;
	int count;
	char buf[BUFSIZ];
	CMD *cmd;

	if ( cmd_id == TGDBLANKLINE )
		return;
	/* one character in argcheck for each argument */
	cmd = cmdbyid( cmd_id );
	for ( s = cmd->argcheck, count = 0; *s; s++, count++ ) {
		if ( count > (argc - 1) )
		{
			tgdputerr( "not enough parameters" );
			return -1;
		}
		switch( *s ) {

		case 's' :	break; /* Anything is a string */

		case 'd' :	t = argv[count];
				if ( *t == '-'  ||  *t == '+' )
					++t;
				for ( ; *t; t++ )
					if ( ! isdigit( *t ) ) {
						sprintf( buf,
						"%s is not an integer",
						argv[count] );
						tgdputerr( buf );
						return -1;
					}
				break;

		case '*' :	goto ENDCHECK; break;

		default :	die( "illegal argcheck" ); break;

		}
	}
	if ( count < argc ) /* should it be <= ? XXX  */
	{
		tgdputerr( "too many parameters" );
		return -1;
	}
	ENDCHECK:
		return 0;
}


/* --
 *  Find the routine to handle the command.
 */
tgdcommand( argc, argv )
int argc;
char *argv[];
{
	char buf[BUFSIZ];
	int cmd_id;
	static int check_syntax_only = 0; /* if 1 do not execute commands */

	cmd_id = which_command( argc, argv );

	/* Metacommands first */
	if ( cmd_id == TGDBLANKLINE  ||  cmd_id == TGDCOMMENT )
		return;
	if  ( cmd_id == TGDCHECKSYNTAX )
		check_syntax_only = 1;
	if ( check_syntax_only ) {
		tgdargcheck( cmd_id, argc, argv ); 
		return;
	}

	/* If only checking syntax, this line is not reached. */

	if ( tgdargcheck( cmd_id, argc, argv ) != 0 )
		die( "parameter syntax mistake" );

	switch ( cmd_id ) {
	case IMAGEARC : imagearc( argc, argv ); break;
	case IMAGEBLUE : imageblue( argc, argv ); break;
	case IMAGEBOUNDSSAFE : imageboundssafe( argc, argv ); break;
	case IMAGECHAR : imagechar( argc, argv ); break;
	case IMAGECHARUP : imagecharup( argc, argv ); break;
	case IMAGECOLORALLOCATE : imagecolorallocate( argc, argv ); break;
	case IMAGECOLORCLOSEST : imagecolorclosest( argc, argv ); break;
	case IMAGECOLORDEALLOCATE : imagecolordeallocate( argc, argv ); break;
	case IMAGECOLOREXACT : imagecolorexact( argc, argv ); break;
	case IMAGECOLORTRANSPARENT : imagecolortransparent( argc, argv ); break;
	case IMAGECOPY : imagecopy( argc, argv ); break;
	case IMAGECOPYRESIZED : imagecopyresized( argc, argv ); break;
	case IMAGECREATE : imagecreate( argc, argv ); break;
	case IMAGECREATEFROMGD : imagecreatefromgd( argc, argv ); break;
	case IMAGECREATEFROMGIF : imagecreatefromgif( argc, argv ); break;
	case IMAGECREATEFROMXBM : imagecreatefromxbm( argc, argv ); break;
	case IMAGECOLORSTOTAL : imagecolorstotal( argc, argv ); break;
	case IMAGEDASHEDLINE : imagedashedline( argc, argv ); break;
	case IMAGEDESTROY : imagedestroy( argc, argv ); break;
	case IMAGEFILL : imagefill( argc, argv ); break;
	case IMAGEFILLTOBORDER : imagefilltoborder( argc, argv ); break;
	case IMAGEFILLEDRECTANGLE : imagefilledrectangle( argc, argv ); break;
	case IMAGEGD : imagegd( argc, argv ); break;
	case IMAGEGETINTERLACED : imagegetinterlaced( argc, argv ); break;
	case IMAGEGETPIXEL : imagegetpixel( argc, argv ); break;
	case IMAGEGETTRANSPARENT : imagegettransparent( argc, argv ); break;
	case IMAGEGIF : imagegif( argc, argv ); break;
	case IMAGEGREEN : imagegreen( argc, argv ); break;
	case IMAGEINTERLACE : imageinterlace( argc, argv ); break;
	case IMAGELINE : imageline( argc, argv ); break;
	case IMAGEFILLEDPOLYGON : imagefilledpolygon( argc, argv ); break;
	case IMAGEPOLYGON : imagepolygon( argc, argv ); break;
	case IMAGERECTANGLE : imagerectangle( argc, argv ); break;
	case IMAGERED : imagered( argc, argv ); break;
	case IMAGESX : imagesx( argc, argv ); break;
	case IMAGESY : imagesy( argc, argv ); break;
	case IMAGESETBRUSH : imagesetbrush( argc, argv ); break;
	case IMAGESETPIXEL : imagesetpixel( argc, argv ); break;
	case IMAGESETSTYLE : imagesetstyle( argc, argv ); break;
	case IMAGESETTILE : imagesettile( argc, argv ); break;
	case IMAGESTRING : imagestring( argc, argv ); break;
	case IMAGESTRINGUP : imagestringup( argc, argv ); break;
	case TGDDEBUG : tgddebug( argc, argv ); break;
	case TGDSTOCKIMAGE : tgdstockimage( argc, argv ); break;
	case TGDSTRINGCENTER : tgdstringcenter( argc, argv ); break;
	default :
		sprintf( buf, "Command %s unrecognized, id = %d.",
			argv[0], cmd_id );
		die ( buf );
		break;
	};
}
/* ---------------------------------------------------------- */

#ifndef TGD_H
#define TGD_H 1

/*
 * tgd.h -- header file for the tgd program
 * SCCS = "@(#)tgd.h	1.1 6/27/95 tgd";
 * Bradley K. Sherman 1995 (bks@netcom.com bks@s27w007.pswfs.gov)
 */
 
#define MAXTIM  64
#define MAXTIMNAME  100
#define MAXCOLOUR 256
#define MAXCOLOURNAME MAXTIMNAME
#define MAXFONTNAME MAXTIMNAME
#define MAXTGDFONT 132
#define TGDMAXLINE 5012
#define TGDMAXTOKENS 100

#define DOUBLE_QUOTE '"'


typedef struct tgdfont {
	gdFontPtr fo;
	char name[MAXFONTNAME];
} TGDFONT; /* Map user supplied font name to a gd font pointer */

typedef struct colour {
	int id;
	char name[MAXCOLOURNAME];
} COLOUR;

typedef struct tim {
	gdImagePtr im;
	char name[MAXTIMNAME];  /* User name for this image */
	COLOUR colour[MAXCOLOUR];  /* Map colour names to gd id's */
} TIM; /* Map user supplied image name to a gd image pointer */

#endif /* TGD_H */ /* End of file tgd.h */

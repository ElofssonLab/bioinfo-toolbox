# Line added by Darrin Stewart to place tgd in /bin
# The change is  under $(PROG)
#
#
#	This is the Makefile for tgd
#
#	tgd is a text interface to the gd library
#
#	Bradley K. Sherman, June 1995
#	Copyright (c) United States of America
#
#	The basic idea here is that you need the
#	gd1.2 library to compile.  This makefile
#	assumes that it is in a subdirectory called
#	gd1.2
#
#	See the Makefile in the gd1.2 subdirectory.
#	You may have to do a ranlib(1) on the gd library
#	on some systems.
#
#	10/11/95 decomposed tgdsamfont.c into separate
#	files to reduce object size.
#	
#	2/9/96 Modified tgdgd.c so that special handle "stdout"
#	results in output of image file to standard output.
#	
#            SCCS = "@(#)Makefile	1.5 2/9/96 tgd"
#
SHELL = /bin/sh
PROG     = tgd
HOME     = .
GDHOME   = $(HOME)/gd1.2
CC	 = cc
CDEFS	 = 
INCLUDE  = -I$(HOME) -I$(GDHOME)
GDLIB    = $(GDHOME)/libgd.a

#
# Some debugging/profiling options
# CFLAGS   = -p  $(CDEFS) $(INCLUDE)
# LDFLAGS 	= -p -n 
# CFLAGS	  = -O $(CDEFS) $(INCLUDE)
# LDFLAGS  = -n 
#
# CFLAGS	 = -g  $(CDEFS) $(INCLUDE)
# LDFLAGS	 = -g -n 
CFLAGS	 = $(CDEFS) $(INCLUDE)
LDFLAGS	 =

LDLIBS	 = -L$(GDHOME) -lgd -lm
MAKEFILE = Makefile
FONT1    = tgdFont10x20.c tgdFont12x24.c tgdFont5x8.c tgdFont6x10.c
FONT2    = tgdFont6x12.c tgdFont6x13.c tgdFont6x13bold.c tgdFont6x9.c
FONT3    = tgdFont7x13.c tgdFont7x13bold.c tgdFont7x14.c tgdFont8x13.c
FONT4    = tgdFont8x13bold.c tgdFont8x16.c tgdFont9x15.c tgdFont9x15bold.c
FONT1O   = tgdFont10x20.o tgdFont12x24.o tgdFont5x8.o tgdFont6x10.o
FONT2O   = tgdFont6x12.o tgdFont6x13.o tgdFont6x13bold.o tgdFont6x9.o
FONT3O   = tgdFont7x13.o tgdFont7x13bold.o tgdFont7x14.o tgdFont8x13.o
FONT4O   = tgdFont8x13bold.o tgdFont8x16.o tgdFont9x15.o tgdFont9x15bold.o
FONTS    = $(FONT1) $(FONT2) $(FONT3) $(FONT4)
FONTSO   = $(FONT1O) $(FONT2O) $(FONT3O) $(FONT4O)
SRC      = tgd.c tgdgd.c tgdcommand.c $(FONTS)
OBJ      = tgd.o tgdgd.o tgdcommand.o $(FONTSO)

$(PROG): $(OBJ) $(GDLIB)
	$(CC) -o $(PROG) $(LDFLAGS) $(OBJ) $(LDLIBS); mv $(PROG) ../../bin

$(OBJ) :  tgd.h
	$(CC)  $(CFLAGS) -c $*.c

$(GDLIB):
	cd $(GDHOME); make

.c.o:
	$(CC) $(CFLAGS) -c $*.c

sources:	$(SRC)

lint:
	lint -x -h -I$(INCLUDE) $(SRC) $(LINTLIBS)

crossref:
	ctags -x $(SRC)




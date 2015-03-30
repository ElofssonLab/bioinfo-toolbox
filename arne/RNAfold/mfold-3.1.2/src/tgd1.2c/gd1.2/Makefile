#Depending on your system, you will need to modify this makefile.

#If you do not have gcc, change the setting for CC, but you must
#use an ANSI standard C compiler (NOT the old SunOS 4.1.3 cc
#compiler; get gcc if you are still using it). 

#If the ar command fails on your system, consult the ar manpage
#for your system. 

# SOLARIS note: delete the ranlib line below!

SHELL = /bin/sh
CC=cc 
AR=ar
CFLAGS=-O
LIBS=-L./ -lgd -lm

all: libgd.a

# ranlib line added by bks for tgd (you might want to remove it)
libgd.a: gd.o gdfontt.o gdfonts.o gdfontmb.o gdfontl.o gdfontg.o \
	gd.h gdfontt.h gdfonts.h gdfontmb.h gdfontl.h gdfontg.h
	rm -f libgd.a
	$(AR) rc libgd.a gd.o gdfontt.o gdfonts.o gdfontmb.o \
		gdfontl.o gdfontg.o
#	ranlib libgd.a


clean:
	rm -f *.o *.a

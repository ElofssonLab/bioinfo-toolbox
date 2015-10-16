# Makefile for matchMSAtopo.cpp
CC         = g++
CFLAGS     = -Wall -O3
LIBS       = -lm 
LIB_PATH   = .
INCLUDE    = .
SRC        = matchMSAtopo.cpp
OBJ        = $(SRC:.cpp=.o) 
EXE        = matchMSAtopo
RM         = /bin/rm -f
CP         = /bin/cp -f
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -L$(LIB_PATH) -o $(EXE)  $(LIBS)
#compile and assemble C++/C source files into object files
# -c flag tells the compiler to create only OBJ files
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $(SRC) 
install:
	$(CP)  $(EXE)  $(BINPATH)
clean:
	$(RM)  $(OBJ)

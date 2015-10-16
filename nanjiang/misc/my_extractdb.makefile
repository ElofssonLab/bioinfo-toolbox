# Makefile for my_extractdb.cpp
CC         = g++
CFLAGS     = -Wall -O3
LIBS       = -lm 
SRC        = my_extractdb.cpp
OBJ        = $(SRC:.cpp=.o) 
EXE        = my_extractdb
RM         = /bin/rm -f
CP         = /bin/cp -f
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ)  -o $(EXE)  $(LIBS)
#compile and assemble C++/C source files into object files
# -c flag tells the compiler to create only OBJ files
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c  $(SRC) 
install:
	$(CP)  $(EXE)  $(BINPATH)
clean:
	$(RM)  $(OBJ)

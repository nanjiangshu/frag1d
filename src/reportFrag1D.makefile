# Makefile for  checkfirst.cpp
CC         = g++
CFLAGS     = -Wall -O3
CFLAGS_DEBUG     = -Wall -g3
LIBS       = -lm -lmyfunc
LIBS_DEBUG       = -lm -lmyfunc_debug
LIB_PATH   = ../lib
INCLUDE    = ../src
SRC        = reportFrag1D.cpp
OBJ        = $(SRC:.cpp=.o) 
OBJ_DEBUG        = $(SRC:.cpp=_debug.o) 
EXE        = reportFrag1D
EXE_DEBUG        = reportFrag1D_debug
RM         = /bin/rm -f
CP         = /bin/cp -f
all: $(EXE)
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -L$(LIB_PATH) -o $(EXE)  $(LIBS)
%.o: %.cpp
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $< -o $@

debug: $(EXE_DEBUG)
$(EXE_DEBUG): $(OBJ_DEBUG)
	$(CC) $(CFLAGS_DEBUG)  $(OBJ_DEBUG) -L$(LIB_PATH) -o $(EXE_DEBUG)  $(LIBS_DEBUG)
%_debug.o: %.cpp
	$(CC) $(CFLAGS_DEBUG) -c -I$(INCLUDE) $< -o $@

install:
	$(CP)  $(EXE)  ../bin/

installdebug:
	$(CP)  $(EXE_DEBUG)  ../bin/

clean:
	$(RM)  $(OBJ) $(EXE)
	$(RM)  $(OBJ_DEBUG) $(EXE_DEBUG)

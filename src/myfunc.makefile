# Makefile for creating the library for Frag1D program, including mypro.cpp,
# mypro.cpp and subsetc.pp 
CC         = g++
CFLAGS     = -c -Wall -O3
CFLAGS_DEBUG     = -c -Wall -g3
FRAG1D   = ../
LDFLAGS    =
MODE       = shared
AR         = ar
ARFLAGS    = rc
LIBS       = -lc
LIB_PATH   = ../lib
INCLUDE    = ../src
SRC        = myfunc.cpp mypro.cpp subset.cpp
OBJ        = $(SRC:.cpp=.o)
OBJ_DEBUG  = $(SRC:.cpp=_debug.o)
SONAME0    = libmyfunc.so
SONAME     = libmyfunc.so.1
SOLIBNAME  = libmyfunc.so.1.0.1
ALIBNAME   = libmyfunc.so.a
SONAME0_DEBUG    = libmyfunc_debug.so
SONAME_DEBUG     = libmyfunc_debug.so.1
SOLIBNAME_DEBUG  = libmyfunc_debug.so.1.0.1
ALIBNAME_DEBUG   = libmyfunc_debug.so.a
CP         = /bin/cp -f
LN         = /bin/ln -sf
RM         = /bin/rm -f

all: $(SOLIBNAME)
$(SOLIBNAME): $(OBJ)
	$(CC) -shared -Wl,-soname,$(SONAME) -o $@ $(OBJ) $(LIBS)

%.o: %.cpp
	$(CC) -fPIC $(CFLAGS) -I$(INCLUDE) $< -o $@

debug: $(SOLIBNAME_DEBUG)
$(SOLIBNAME_DEBUG): $(OBJ_DEBUG)
	$(CC) -shared -Wl,-soname,$(SONAME) -o $@ $(OBJ_DEBUG) $(LIBS)

%_debug.o: %.cpp
	$(CC) -fPIC $(CFLAGS_DEBUG) -I$(INCLUDE) $< -o $@


install:
	$(CP) $(SOLIBNAME) $(LIB_PATH)
	cd ${LIB_PATH} ; $(LN) $(SOLIBNAME) ${SONAME}; $(LN) $(SONAME) $(SONAME0)

installdebug:
	$(CP) $(SOLIBNAME_DEBUG) $(LIB_PATH)
	cd ${LIB_PATH} ; $(LN) $(SOLIBNAME_DEBUG) ${SONAME_DEBUG}; $(LN) $(SONAME_DEBUG) $(SONAME0_DEBUG)

clean:
	$(RM)  $(OBJ)  $(SOLIBNAME)
	$(RM)  $(OBJ_DEBUG)  $(SOLIBNAME_DEBUG)

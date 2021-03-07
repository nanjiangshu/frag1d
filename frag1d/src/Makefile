# Makefile for compiling Frag1D
export FRAG1D := ${PWD}/..
export LD_LIBRARY_PATH := ../lib
CP = /bin/cp -f 

all:
	make -f myfunc.makefile
	make install -f myfunc.makefile
	make -f rootname.makefile
	make -f pssm2Qij.makefile
	make -f pssm2modm.makefile
	make -f checkfirst.makefile
	make -f checkresult.makefile
	make -f search_new.makefile
	make -f build_hsrfrag.makefile
	make -f reportFrag1D.makefile
debug:
	make debug -f myfunc.makefile
	make installdebug -f myfunc.makefile
	make debug -f rootname.makefile
	make debug -f pssm2Qij.makefile
	make debug -f pssm2modm.makefile
	make debug -f checkfirst.makefile
	make debug -f checkresult.makefile
	make debug -f search_new.makefile
	make debug -f build_hsrfrag.makefile
	make debug -f reportFrag1D.makefile
install:
	make install -f rootname.makefile
	make install -f pssm2Qij.makefile
	make install -f pssm2modm.makefile
	make install -f checkfirst.makefile
	make install -f checkresult.makefile
	make install -f search_new.makefile
	make install -f build_hsrfrag.makefile
	make install -f reportFrag1D.makefile
clean:
	make clean  -f myfunc.makefile
	make clean  -f rootname.makefile
	make clean  -f pssm2Qij.makefile
	make clean  -f pssm2modm.makefile
	make clean  -f checkfirst.makefile
	make clean  -f checkresult.makefile
	make clean  -f search_new.makefile
	make clean  -f build_hsrfrag.makefile
	make clean  -f reportFrag1D.makefile

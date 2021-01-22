#ROOT includes + libraries
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)


#Debugging info + allWarnings + level-3 optimalisation
CFLAGS=-g -Wall -O3

all: example-6 example-7 

# Rule to compile examples
example-%: example-%.cc ranlxd.o
	g++ $(CFLAGS) $(ROOTCFLAGS)   $^  $(ROOTLIBS)  -o $@ 


# Rule to compile Ranlux
ranlxd.o: ranlxd.cc
	g++ $(CFLAGS) -DSSE -c  $^ -o $@  

clean:
	rm -rf *.o example-? example*.dSYM



#LHAPDF includes + libraries
isLHA := $(shell which lhapdf-config > /dev/null; echo $$?)
ifeq ($(isLHA),0)
LHALIBS = $(shell lhapdf-config --libdir)
LHAINC  = $(shell lhapdf-config --incdir)
else
LHAdir=/home/terascale/lhapdf/local
LHALIBS = $(LHAdir)/lib
LHAINC  = $(LHAdir)/include
endif

isTMDlib := $(shell which ~/jung/cvs/TMDlib/TMDlib2/local/bin/TMDlib-config > /dev/null; echo $$?)
ifeq ($(isTMDlib),0)
TMDLIB = $(shell ~/jung/cvs/TMDlib/TMDlib2/local/bin/TMDlib-config --libdir)
TMDLIBINC=$(shell ~/jung/cvs/TMDlib/TMDlib2/local/bin/TMDlib-config --incdir)/tmdlib
GSLlib=$(shell gsl-config --libs)
else
TMDLIB = /home/terascale/TMDlib/local/lib -lTMDlib
TMDLIBINC=/home/terascale/TMDlib/local//include/tmdlib
GSLlib=$(shell gsl-config --libs)
endif


example-lhapdf: example-lhapdf.cc ranlxd.o
	g++  $(CFLAGS) $(ROOTCFLAGS) -I$(LHAINC)  $^  $(ROOTLIBS) -L$(LHALIBS) -lLHAPDF  -o $@ 

example-tmdlib: example-tmdlib.cc ranlxd.o
	g++  $(CFLAGS) $(ROOTCFLAGS) -I$(LHAINC) -I$(TMDLIBINC) $^  $(ROOTLIBS) -L$(TMDLIB) -lTMDlib -L$(LHALIBS) -lLHAPDF $(GSLlib) -o $@

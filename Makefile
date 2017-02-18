###############################################################################

#include $(ROOTSYS)/etc/Makefile.arch

###############################################################################
OBJ = norm.exe
MainFile = snMain.cc 
###############################################################################

SourceFile := $(wildcard src/*.cc) 
IncludeFile := $(wildcard include/*.hh) 

###############################################################################

ROOTCFLAGS  = $(shell root-config --cflags)
#ROOTLIBS    = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)

GXX = g++ -std=c++0x 
DIR_INC = -Iinclude 
#-I$(TARTSYS)/include
CFLAGS = -Wall -O2 $(DIR_INC) -lTMVA -lMathMore 
#-fopenmp 
#-L$(TARTSYS)/lib
#-lXMLParser -lanacore -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
###############################################################################

all:$(OBJ)
$(OBJ): $(MainFile) $(SourceFile)
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) $(ROOTGLIBS) -o $@ $(MainFile) $(SourceFile)
	@echo "================================"
	@echo "Compile $@ done !"
	@echo "================================"

clean:
	rm -f *.o *.d *Dict.* *.root *.txt $(OBJ)
	
###############################################################################


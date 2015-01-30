#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

#ROOFITINCLUDE = 
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif

CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS)
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lEG -lGenVector
COREDIR = CORE
SSCOREDIR = SSCORE

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lGenVector

#DIR = /Users/cerati/SSAnalysis/SSAnalysis/
DIR = ./

CORESOURCES=$(DIR)/$(COREDIR)/CMS3.cc \
 $(DIR)/$(COREDIR)/ElectronSelections.cc \
 $(DIR)/$(COREDIR)/MuonSelections.cc \
 $(DIR)/$(COREDIR)/JetSelections.cc \
 $(DIR)/$(COREDIR)/MCSelections.cc \
 $(DIR)/$(COREDIR)/MT2/MT2Utility.cc \
 $(DIR)/$(COREDIR)/MT2/MT2.cc \
 $(DIR)/$(COREDIR)/SSSelections.cc #\
#./$(COREDIR)/conversionTools.cc \
#./$(COREDIR)/electronSelectionsParameters.cc \
#./$(COREDIR)/eventSelections.cc \
#./$(COREDIR)/jetSelections.cc \
#./$(COREDIR)/mcSelections.cc \
#./$(COREDIR)/metSelections.cc \
#./$(COREDIR)/MITConversionUtilities.cc \
#./$(COREDIR)/trackSelections.cc \
#./$(COREDIR)/triggerUtils.cc \
#./$(COREDIR)/utilities.cc 
COREOBJECTS=$(CORESOURCES:.cc=.o)
CORELIB=libCORE.so

SSCORESOURCES=$(DIR)/$(SSCOREDIR)/selections.cc 
SSCOREOBJECTS=$(SSCORESOURCES:.cc=.o)
SSCORELIB=libSSCORE.so

SOURCES = $(wildcard $(DIR)/*.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIB = liblooper.so

DICT = LinkDef_out.o

LIBS = $(LIB)

EXE = main.exe

#
# how to make it
#

$(CORELIB): $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

$(SSCORELIB): $(SSCOREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) $(LINKERFLAGS) -shared $(SSCOREOBJECTS) -o $@"; \
	$(LINKER) $(LINKERFLAGS) -shared $(SSCOREOBJECTS) -o $@

$(LIB):	$(DICT) $(OBJECTS) $(COREOBJECTS) $(SSCOREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) $(CFLAGS) $(LINKERFLAGS) -shared $(OBJECTS) $(DICT) -o $@"; \
	$(LINKER) -shared -o $@ $(OBJECTS) $(COREOBJECTS) $(SSCOREOBJECTS) $(DICT) $(LINKERFLAGS)

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(DICTINCLUDE)  LinkDef.h; \
	cat LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

%.exe:  $(LIBS) 
	$(QUIET) echo "Building $@"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) ${@:.exe=.cc} 

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o:    %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

loopclean:
	rm -f \
	*_out.*	 \
	*.o \
	*.*~ \
	$(LIB) \
	main.exe

clean: loopclean
	rm -f \
	$(CORELIB) \
	./$(COREDIR)/*.o \
	./$(COREDIR)/MT2/*.o \
	$(SSCORELIB) \
	./$(SSCOREDIR)/*.o

endif

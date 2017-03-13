#OS_DEFS SOEXE,SOFLAGS,DEFS
PACKAGE=$(shell basename `pwd`)
UNAME= $(shell uname)
ifeq ($(UNAME), Linux)
SOFLAGS     = -shared
DEFS        = -DOS_LINUX
LIB_EXT    := .so
endif
ifeq ($(UNAME), Darwin)
SOFLAGS     = -dynamiclib
DEFS        = -DOS_DARWIN
LIB_EXT    := .dylib
endif
#AR           = ar
#LD_SHARED    = g++

CXX	     = g++
CC           = gcc
ROOTCINT     = $(ROOTSYS)/bin/rootcint
DEFS		+= -DCERNLIB_TYPE -DCERNLIB_LINUX -DCERNLIB_UNIX -DHEP_SHORT_NAMES
GENFITDIR    =$(GENFIT_ROOT_DIR)
#/Users/jwlee/local/util/genfit/GenFit_INSTALL
BOOSTDIR     = $(BOOST_ROOT_DIR)
#/Users/jwlee/local/util/boost/boost_1_61_0_compile

CXX_SOURCES  = $(wildcard src/*.cc)
CXX_DICT     = $(wildcard include/*LinkDef.h)
CXX_DICTSRC  = $(CXX_DICT:include/%LinkDef.h=src/%Dict.cxx)
CXX_OBJECTS  = $(filter-out %Dict.o,$(CXX_SOURCES:src/%.cc=tmp/%.o))
CXX_DICTOBJ  = $(CXX_DICTSRC:src/%.cxx=tmp/*.o)

CXXFLAGS     = -fPIC -Wall
CXXFLAGS    += $(shell $(ROOTSYS)/bin/root-config --cflags)
CXXFLAGS    += -g -fmessage-length=0 -fpermissive #-O2

LIBS        += $(shell $(ROOTSYS)/bin/root-config --glibs) -lEG -lGeom
LIBS        += -L/usr/local/Cellar/libusb/1.0.20/lib -lusb-1.0
LIBS        += -L$(GENFITDIR)/lib -lgenfit2
LIBS        += -L$(BOOSTDIR)/lib
LIBS        += -L$(BOOSTDIR)/libs -lboost_atomic -lboost_filesystem -lboost_math_c99 -lboost_prg_exec_monitor -lboost_signals -lboost_wave -lboost_chrono -lboost_graph -lboost_math_c99f -lboost_program_options -lboost_system -lboost_wserialization -lboost_container -lboost_math_c99l -lboost_thread -lboost_context -lboost_locale -lboost_math_tr1 -lboost_random -lboost_timer -lboost_coroutine -lboost_log -lboost_math_tr1f -lboost_regex -lboost_type_erasure -lboost_date_time -lboost_log_setup -lboost_math_tr1l -lboost_serialization -lboost_unit_test_framework #-lboost_python
LIBS        += -L$(E42_TOP_DIR)/lib/so -lGsimData


INCLUDEDIR  += .  ./include
INCLUDEDIR  += $(ROOTSYS)/include
INCLUDEDIR  += $(LIBUSB_INC)
INCLUDEDIR  += $(GENFITDIR)/include
INCLUDEDIR  += $(BOOSTDIR)/include
INCLUDEDIR  += $(BOOSTDIR)
INCLUDEDIR  += $(E42_TOP_DIR)/include

ifeq ($(UNAME), Linux)
INCLUDEDIR  += /usr/include/libusb-1.0
endif
ifeq ($(UNAME), Darwin)
INCLUDEDIR  += /usr/local/Cellar/libusb/1.0.20/include/libusb-1.0
endif

INCLUDES     = $(INCLUDEDIR:%=-I%)

COMPILE_CC  := $(CC)  -c $(DEFS) $(INCLUDES) $(CFLAGS)
COMPILE_CXX := $(CXX) -c $(DEFS) $(INCLUDES) $(CXXFLAGS)
DEPEND_CC   := $(CC) -MM $(DEFS) $(INCLUDES) $(filter-out -fPIC, $(CFLAGS))
DEPEND_CXX  := $(CC) -MM $(DEFS) $(INCLUDES) $(filter-out -fPIC, $(CXXFLAGS))
LINK_CC     := $(CC)
LINK_CXX    := $(CXX)

TARGET_SRC   = $(wildcard *.cc)
TARGET_BIN   = $(TARGET_SRC:%.cc=bin/%)
#TARGET_OBJS  = $(TARGET_SRC:%.cc=tmp/%.o)

all: rootDict makelib target

target: $(TARGET_BIN)
	@ basename `pwd`
makelib: $(CXX_OBJECTS)
	@ echo $(ROOTSYS)
	@ echo $(CXX_OBJECTS)
	@[ -d ./lib ] || mkdir -p ./lib
	$(CXX) $(SOFLAGS) -o lib/lib$(PACKAGE)$(LIB_EXT) $(CXX_OBJECTS) $(CXX_DICTSRC) $(CXXFLAGS) $(INCLUDES) $(LIBS)
	@ echo "make library"
rootDict:
ifneq "$(shell ls ./include/*LinkDef.h 2>/dev/null)" ""
	@ echo making rootdicts
	@ for i in `ls ./include/*LinkDef.h | sed -n "s/\(.*\)LinkDef.h/\1.h/p"`;do \
	j=`basename $$i .h`; \
	$(ROOTCINT) -f src/$${j}Dict.cxx -c $(INCLUDES) ./include/$${j}.h ./include/$${j}LinkDef.h; \
	echo $$i ./include/$${j}LinkDef.h; \
	done;
	@ echo End making rootdicts
endif

bin/%: tmp/%.o $(CXX_OBJECTS)
	@ echo making $(@F)
	@ [ -d ./bin ] || mkdir -p ./bin
	@ echo "##########################################################"
	@ echo COMPILE $(@F)
	@ echo "##########################################################"
	$(CXX) -o $@ $^ $(filter-out -fPIC,$(CXXFLAGS)) $(CXX_DICTSRC) $(INCLUDES) $(LIBS)
	@ echo "##########################################################"

tmp/%.o: src/%.c
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making $(@F) with $(<F) and $(^F)
	@ $(COMPILE_CC) -c $< $(CXXFLAGS) -o $@

tmp/%.d: src/%.c
	@ [-d ./tmp ] || mkdir -p ./tmp
	@ echo making dependencies of $(<F)
	@ $(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*\.o[ :]*/tmp\/$(@F) tmp\/&/g" > $@'

tmp/%.o: src/%.cc
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making $(@F) with $(<F) and $(^F)
	@ $(COMPILE_CXX) $<  $(INCLUDES) -o $@

tmp/%.d: src/%.cc
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making dependencies of $(<F)
	@ $(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*\.o[ :]*/tmp\/$(@F) tmp\/&/g" > $@'


tmp/%.o: %.cc
	@ echo $(INCLUDEDIR)
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making $(@F)
	$(COMPILE_CXX) $< -o $@

tmp/%.d: %.cc
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making dependencies of $(@F)
	@ $(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*\.o[ :]*/tmp\/$(@F) tmp\/&/g" > $@'




.PHONY: clean
clean:
	rm -f $(TARGET_BIN)
	rm -f ./lib/*
	rm -f ./tmp/*.d
	rm -f ./src/*Dict.*
	rm -f ./tmp/*.o
	rm -f $(TARGET_OBJS)

# Makefile for divstats
#
#   make              build ./divstats
#   make check        build, then run the golden-file regression suite
#   make clean        remove objects and the binary
#   make install      install to $(PREFIX)/bin   (default /usr/local)
#   make info         print the detected platform and the flags in use
#
# The platform is DETECTED from uname; nothing here needs editing to build.
# Every variable below is overridable on the command line, e.g.
#
#   make CXX=clang++ OPT="-O2 -g"
#   make L_PATH=../lib/linux          # unusual layout or cross build
#   make install PREFIX=$HOME/.local
#
# Previously the three platform blocks were selected by commenting and
# uncommenting, with a second near-duplicate file (Makefile.lin) carrying the
# same content with different lines commented -- so a change to the build had
# to be made twice, and the default block was the wrong one for most readers.

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# `CXX ?= g++` does nothing: make defines CXX = c++ via its builtin rules, so
# the variable is already set and ?= never fires. Override make's DEFAULT
# while still respecting an explicit command-line or environment value.
ifeq ($(origin CXX),default)
  CXX := g++
endif
OPT      ?= -O3

# An explicit dialect floor. With no -std= flag the dialect was whatever the
# compiler happened to default to, which drifts between compiler versions and
# silently changes which language features compile.
STD      ?= -std=c++11

# -m64 -mmmx -msse -msse2 are x86 codegen flags. They do not break an ARM
# build -- clang accepts and ignores them -- but they emit
# "argument unused during compilation" on every translation unit, which is
# noise that trains the reader to skim warnings. Apply them only where they
# mean something.
ifeq ($(filter x86_64 amd64,$(UNAME_M)),)
  ARCHFLAGS ?=
else
  ARCHFLAGS ?= -m64 -mmmx -msse -msse2
endif

I_PATH   ?= -I../include
LINK_OPTS2 ?= -lz -pthread

# Where the vendored static libhts.a for this platform lives.
ifeq ($(UNAME_S),Darwin)
  ifeq ($(UNAME_M),arm64)
    L_PATH ?= ../lib/macos-arm
  else
    L_PATH ?= ../lib/osx
  endif
else
  ifeq ($(UNAME_S),Linux)
    L_PATH ?= ../lib/linux
  else
    # MinGW reports MINGW64_NT-10.0 and friends
    L_PATH ?= ../lib/win32
    I_PATH ?= -I../include -I../include/win32
    LINK_OPTS2 ?= ../lib/win32/libpthreadGC2.a ../lib/win32/libz.a
    OPT ?= -DPTW32_STATIC_LIB -O3 -static-libgcc -static-libstdc++
  endif
endif

CXXFLAGS ?= $(OPT) $(STD) $(ARCHFLAGS)

PREFIX  ?= /usr/local
INSTALL ?= install

#htslib supplies the VCF/BCF reader. Configured without bzip2, lzma, libcurl
#and libdeflate (those are CRAM features), so it needs nothing beyond the -lz
#already linked -- see lib/build_htslib.sh.
#
#GSL is gone: gsl_combination was its only use, in the pi_k2 tie enumeration
#that the closed form replaces.
LINK_OPTS = $(L_PATH)/libhts.a

#Every object depends on every header. Over-broad but always correct: before
#this, the .o rules had no header dependencies at all, so changing a header
#left stale objects linked against the old declarations. That is not a
#cosmetic issue -- changing segsites() from int to double while
#divstats-wintools.o was stale produced a binary that read the return value
#from the integer register and reported S = 5.66e+07.
HEADERS = divstats-cli.h divstats-data.h divstats-hts.h divstats-winstats.h divstats-wintools.h \
          param_t.h binom.h gzstream.h

OBJS = divstats-main.o divstats-data.o divstats-hts.o binom.o param_t.o \
       gzstream.o divstats-winstats.o divstats-wintools.o

all : divstats

#A missing vendored archive used to surface as a wall of undefined symbols
#from the linker. Say what is missing and how to produce it instead.
$(L_PATH)/libhts.a :
	@echo "ERROR: $(L_PATH)/libhts.a is missing."
	@echo "       divstats vendors a static htslib per platform. Build one with:"
	@echo "           ../lib/build_htslib.sh"
	@echo "       Detected platform: $(UNAME_S) $(UNAME_M)"
	@exit 1

divstats : $(OBJS) $(L_PATH)/libhts.a
	$(CXX) $(CXXFLAGS) -o divstats $(OBJS) $(LINK_OPTS2) $(LINK_OPTS)

divstats-main.o : divstats-main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c divstats-main.cpp $(I_PATH)

divstats-data.o : divstats-data.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c divstats-data.cpp $(I_PATH)

divstats-hts.o : divstats-hts.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c divstats-hts.cpp $(I_PATH)

divstats-winstats.o : divstats-winstats.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c divstats-winstats.cpp $(I_PATH)

divstats-wintools.o : divstats-wintools.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c divstats-wintools.cpp $(I_PATH)

binom.o : binom.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c binom.cpp

param_t.o : param_t.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c param_t.cpp

gzstream.o : gzstream.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c gzstream.cpp $(I_PATH)

#Golden-file regression suite. Compares this build's output against the
#tables in tests/expected/ column by column. See tests/run_tests.sh.
check : divstats
	../tests/run_tests.sh

clean :
	rm -f $(OBJS) divstats

install : divstats
	$(INSTALL) -d $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -m 755 divstats $(DESTDIR)$(PREFIX)/bin/divstats

#Printed on request, and worth pasting into a bug report.
info :
	@echo "platform : $(UNAME_S) $(UNAME_M)"
	@echo "CXX      : $(CXX)"
	@echo "CXXFLAGS : $(CXXFLAGS)"
	@echo "L_PATH   : $(L_PATH)"
	@echo "libhts.a : $(if $(wildcard $(L_PATH)/libhts.a),present,MISSING -- run ../lib/build_htslib.sh)"

.PHONY : all check clean install info

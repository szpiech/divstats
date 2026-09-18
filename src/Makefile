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

# Windows, built with the MinGW-w64 toolchain under MSYS2. uname -s reports
# MINGW64_NT-10.0-<build> from a MINGW64 shell, MSYS_NT-10.0-<build> from the
# MSYS shell, and CYGWIN_NT-* under Cygwin, so match on the prefix rather than
# enumerating build numbers. This is a real check now: it used to be the
# catch-all else, so every unrecognised platform -- FreeBSD, Solaris, a
# cross-compile -- was silently handed Windows library paths.
WINDOWS := $(filter MINGW% MSYS% CYGWIN%,$(UNAME_S))

# Link-time executable suffix, empty everywhere but Windows so the target
# below stays plain `divstats` on Unix. Not cosmetic: MinGW's g++ appends .exe
# whatever you pass to -o, so without this the file named by the target is
# never created, make considers the goal permanently out of date, and every
# invocation relinks.
ifneq ($(WINDOWS),)
  EXE := .exe
else
  EXE :=
endif

# Where the vendored static libhts.a for this platform lives.
ifneq ($(WINDOWS),)
  WINFLAGS ?= -static-libgcc -static-libstdc++
  L_PATH ?= ../lib/mingw64
  # MSYS2's mingw-w64 toolchain supplies zlib and winpthreads, so the two
  # 32-bit archives lib/win32 used to carry are neither needed nor usable from
  # a 64-bit build.
  LINK_OPTS2 ?= -lz -pthread
  # Static libgcc/libstdc++ so the .exe runs outside the MSYS2 shell, where
  # those DLLs are not on PATH. -DPTW32_STATIC_LIB went with pthreads-win32;
  # winpthreads needs no such define.
  #
  # Its own variable, not OPT: `OPT ?= -O3` is set above this block, so an
  # `OPT ?=` here can never fire -- which is why the old win32 branch's
  # -static-libgcc and -DPTW32_STATIC_LIB were dead the whole time. It also
  # keeps the flags when a caller overrides OPT, as the sanitizer and warning
  # CI jobs do.
else
  ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),arm64)
      L_PATH ?= ../lib/macos-arm
    else
      L_PATH ?= ../lib/osx
    endif
  else
    ifeq ($(UNAME_S),Linux)
      # lib/linux holds the x86_64 archive; arm64 needs its own, exactly as
      # Darwin does. An archive built for one architecture cannot link against
      # objects compiled for the other, and the failure is a wall of undefined
      # symbols rather than anything that names the cause.
      #
      # uname -m reports aarch64 on most arm64 distributions and arm64 on a
      # few, so match both.
      ifeq ($(filter aarch64 arm64,$(UNAME_M)),)
        L_PATH ?= ../lib/linux
      else
        L_PATH ?= ../lib/linux-arm
      endif
    else
      # Deliberately a path that cannot exist. `make info` still answers, and
      # the missing-archive rule then names the platform and says to pass
      # L_PATH -- which beats compiling against another platform's headers.
      L_PATH ?= ../lib/unknown-$(UNAME_S)-$(UNAME_M)
    endif
  endif
endif

CXXFLAGS ?= $(OPT) $(STD) $(ARCHFLAGS) $(WINFLAGS)

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

all : divstats$(EXE)

#A missing vendored archive used to surface as a wall of undefined symbols
#from the linker. Say what is missing and how to produce it instead.
$(L_PATH)/libhts.a :
	@echo "ERROR: $(L_PATH)/libhts.a is missing."
	@echo "       divstats vendors a static htslib per platform. Build one with:"
	@echo "           ../lib/build_htslib.sh"
	@echo "       Detected platform: $(UNAME_S) $(UNAME_M)"
	@echo "       If divstats does not know that platform, pass the path"
	@echo "       explicitly: make L_PATH=../lib/<dir>"
	@exit 1

divstats$(EXE) : $(OBJS) $(L_PATH)/libhts.a
	$(CXX) $(CXXFLAGS) -o divstats$(EXE) $(OBJS) $(LINK_OPTS2) $(LINK_OPTS)

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
check : divstats$(EXE)
	../tests/run_tests.sh

clean :
	rm -f $(OBJS) divstats divstats.exe

install : divstats$(EXE)
	$(INSTALL) -d $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -m 755 divstats$(EXE) $(DESTDIR)$(PREFIX)/bin/divstats$(EXE)

#Printed on request, and worth pasting into a bug report.
info :
	@echo "platform : $(UNAME_S) $(UNAME_M)"
	@echo "CXX      : $(CXX)"
	@echo "CXXFLAGS : $(CXXFLAGS)"
	@echo "L_PATH   : $(L_PATH)"
	@echo "binary   : divstats$(EXE)"
	@echo "libhts.a : $(if $(wildcard $(L_PATH)/libhts.a),present,MISSING -- run ../lib/build_htslib.sh)"

.PHONY : all check clean install info

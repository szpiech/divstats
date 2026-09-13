# Makefile for divstats

#debugging
#G++FLAG = -g -Wall -Werror
#For dynamic linking
#LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas

#for osx systems
CC = g++
G++FLAG = -O3 -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz -lgsl -lgslcblas
LINK_OPTS2 = -lz -pthread
#-lgsl -lgslcblas
I_PATH = -I../include
#L_PATH = ../lib/osx
L_PATH = ../lib/macos-arm

#for linux systems
#CC = g++
#G++FLAG = -O3 -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz
#I_PATH = -I../include
#L_PATH = ../lib/linux

#for windows, using MinGW build environment
#CC = g++.exe
#G++FLAG = -DPTW32_STATIC_LIB -O3 -static-libgcc -static-libstdc++
#LINK_OPTS2 = ../lib/win32/libpthreadGC2.a ../lib/win32/libz.a 
#I_PATH = -I../include -I../include/win32 
#L_PATH = ../lib/win32

#For static linking of norm program to gsl libs
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

all : divstats

divstats : divstats-main.o divstats-data.o divstats-hts.o binom.o param_t.o gzstream.o divstats-winstats.o divstats-wintools.o
	$(CC) $(G++FLAG) -o divstats divstats-main.o divstats-data.o divstats-hts.o binom.o param_t.o gzstream.o divstats-winstats.o divstats-wintools.o $(LINK_OPTS2) $(LINK_OPTS)

divstats-main.o : divstats-main.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c divstats-main.cpp $(I_PATH)

divstats-data.o : divstats-data.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c divstats-data.cpp $(I_PATH)

divstats-hts.o : divstats-hts.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c divstats-hts.cpp $(I_PATH)

divstats-winstats.o : divstats-winstats.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c divstats-winstats.cpp $(I_PATH)

divstats-wintools.o : divstats-wintools.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c divstats-wintools.cpp $(I_PATH)

binom.o : binom.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c binom.cpp

param_t.o : param_t.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.cpp $(HEADERS)
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH)

#Golden-file regression suite. Compares this build's output against the
#tables in tests/expected/ column by column. See tests/run_tests.sh.
check : divstats
	../tests/run_tests.sh

clean :
	rm *.o

# Makefile for divstats

#debugging
#G++FLAG = -g -Wall -Werror
#For dynamic linking
#LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas

#for osx systems
CC = g++
G++FLAG = -fast -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz
LINK_OPTS2 = -lz
I_PATH = -I../include
L_PATH = ../lib/osx

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
#LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a

all : divstats

divstats : divstats-main.o divstats-data.o binom.o param_t.o gzstream.o hamming_t.o divstats-winstats.o
	$(CC) $(G++FLAG) -o divstats divstats-main.o divstats-data.o binom.o param_t.o gzstream.o hamming_t.o divstats-winstats.o $(LINK_OPTS2)

divstats-main.o : divstats-main.cpp
	$(CC) $(G++FLAG) -c divstats-main.cpp $(I_PATH)

divstats-data.o : divstats-data.cpp
	$(CC) $(G++FLAG) -c divstats-data.cpp $(I_PATH)

divstats-winstats.o : divstats-winstats.cpp
	$(CC) $(G++FLAG) -c divstats-winstats.cpp $(I_PATH)

binom.o : binom.cpp
	$(CC) $(G++FLAG) -c binom.cpp

param_t.o : param_t.cpp
	$(CC) $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.cpp
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH)

hamming_t.o : hamming_t.cpp
	$(CC) $(G++FLAG) -c hamming_t.cpp $(I_PATH)

clean :
	rm *.o

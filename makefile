CC=gcc
GROMACS=$(HOME)/gromacs-src/
DEBUG=-ggdb
NO_VALGRIND=-g -fno-inline
NO_PROFILE=-pg -fno-omit-frame-pointer
NO_OMP=-DUSE_OMP
TP=-DTIME_PROFILE
#
#OPT=-O3 -Wall -pedantic -std=c99  -Werror -Wno-unused -Wuninitialized  -ffast-math -fopenmp
OPT=-O -fno-inline -gfull
all: g_density 

clean: 
	rm *.o g_density 
g_density: g_density.c gmx_density.o makefile
	$(CC) $(PROFILE) $(TP) -DNO_ALPHA $(VALGRIND) $(DEBUG) $(OPT)  -L$(HOME)/lib/ -DHAVE_CONFIG_H -I. $(ALPHA_INCLUDE) -I$(GROMACS)/include -I$(GROMACS)/src -I/usr/include/libxml2     -I./include  -o g_density gmx_density.o $(ALPHA_LIB) -ldl -lz -lgmx -lmd  g_density.c 
gmx_density.o: gmx_density.c makefile
	$(CC) $(PROFILE) $(TP) $(OMP) -DNO_ALPHA  $(VALGRIND) $(DEBUG)  $(OPT) -DHAVE_CONFIG_H  -I$(GROMACS)/include -I$(GROMACS)/src -I/usr/include/libxml2     -I./include -c -o gmx_density.o gmx_density.c

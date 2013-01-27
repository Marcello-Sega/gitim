GROMACS=$(HOME)/gromacs-src/
CC=gcc
NO_DEBUG=-ggdb
NO_VALGRIND=-g -fno-inline
NO_PROFILE=-pg -fno-omit-frame-pointer
NO_OMP=-DUSE_OMP
TP=-DTIME_PROFILE
ALPHA_INCLUDE=-I./qhull/
ALPHA_LIB=-L./qhull -lqhull
ALPHA_LIB=-L$(HOME)/Dropbox/Documents/ITIM/qhull_src -lqhull
GMXLIBS=-L$(GROMACS)/src/mdlib/.libs/ -L$(GROMACS)/src/gmxlib/.libs/
#
OPT=-O3 -Wall -pedantic -std=c99  -Werror -Wno-unused -Wuninitialized  -ffast-math -fopenmp
NO_OPT=-O -fno-inline -gfull


all: gitim

distclean:
	(cd qhull; make clean) ; make clean

clean: 
	rm -f *.o gitim 
gitim:  makefile qhull/libqhull.a gitim.c
	$(CC)  $(PROFILE) $(TP) -DALPHA $(VALGRIND) $(DEBUG)  $(OPT)  -DHAVE_CONFIG_H -I. $(ALPHA_INCLUDE) -I$(GROMACS)/include -I$(GROMACS)/src -I/usr/include/libxml2  -I./include -o gitim    $(ALPHA_LIB) $(GMXLIBS) -ldl -lz -lgmx -lm  gitim.c

qhull/libqhull.a: 
	(cd qhull ; make) 

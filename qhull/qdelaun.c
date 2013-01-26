/*<html><pre>  -<a                             href="qh-qhull.htm"
  >-------------------------------</a><a name="TOP">-</a>

   qdelaun.c
     compute Delaunay triangulations and furthest-point Delaunay
     triangulations using qhull

   see unix.c for full interface

   copyright (c) 1993-2010, The Geometry Center
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "libqhull.h"
#include "mem.h"
#include "qset.h"
char hidden_options[]=" d n v H U Qb QB Qc Qf Qg Qi Qm Qr QR Qv Qx TR E V FC Fi Fo Ft Fp FV Q0 Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 ";

/*-<a                             href="qh-qhull.htm#TOC"
  >-------------------------------</a><a name="prompt">-</a>

  qh_prompt
    long prompt for qhull

  notes:
    restricted version of libqhull.c

  see:
    concise prompt below
*/


/* for opts, don't assign 'e' or 'E' to a flag (already used for exponent) */

/*-<a                             href="qh-qhull.htm#TOC"
  >-------------------------------</a><a name="prompt2">-</a>

  qh_prompt2
    synopsis for qhull
*/
/* for opts, don't assign 'e' or 'E' to a flag (already used for exponent) */

/*-<a                             href="qh-qhull.htm#TOC"
  >-------------------------------</a><a name="prompt3">-</a>

  qh_prompt3
    concise prompt for qhull
*/

/*-<a                             href="qh-qhull.htm#TOC"
  >-------------------------------</a><a name="main">-</a>

  main( argc, argv )
    processes the command line, calls qhull() to do the work, and exits

  design:
    initializes data structures
    reads points
    finishes initialization
    computes convex hull and other structures
    checks the result
    writes the output
    frees memory
*/
int main(int argc, char *argv[]) {
  int curlong, totlong; /* used !qh_NOmem */
  int exitcode, numpoints, dim;
  coordT *points;
  boolT ismalloc;


 qh_option("delaunay  incidence Qbbound-last", NULL, NULL);
 qh DELAUNAY= True;     /* 'd'   */
 qh SCALElast= True;    /* 'Qbb' */
 qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */
/**/
 qh TRIangulate = True; /* 'Qt' */
 qh ATinfinity = True ; /* 'Qz' */

 qh_checkflags(qh qhull_command, hidden_options);
 qh_initflags(qh qhull_command);
 points=(coordT*) malloc(3*sizeof(coordT)*4); /* qh_readpoints(&numpoints, &dim, &ismalloc); */
{int n=0;
 points[3*n+0]= -0.5 ; points[3*n+1]=  -0.5  ; points[3*n+2]=   -0.5;n++; 
 points[3*n+0]= -0.5 ; points[3*n+1]=  -0.5  ; points[3*n+2]=    0.5;n++; 
 points[3*n+0]= -0.5 ; points[3*n+1]=   0.5  ; points[3*n+2]=   -0.5;n++; 
 points[3*n+0]= -0.5 ; points[3*n+1]=   0.5  ; points[3*n+2]=    0.5;n++; 
 points[3*n+0]=  0.5 ; points[3*n+1]=  -0.5  ; points[3*n+2]=   -0.5;n++; 
 points[3*n+0]=  0.5 ; points[3*n+1]=  -0.5  ; points[3*n+2]=    0.5;n++; 
 points[3*n+0]=  0.5 ; points[3*n+1]=   0.5  ; points[3*n+2]=   -0.5;n++; 
 points[3*n+0]=  0.5 ; points[3*n+1]=   0.5  ; points[3*n+2]=    0.5;n++; 
}
 dim=3;numpoints=8;ismalloc=True;
 qh ferr = stderr; /*qh_FILEstderr; */
 qh_init_B(points, numpoints, dim, ismalloc);
 exit(printf("here\n"));

 qh_qhull();
 qh_check_output();
 qh_produce_output();
 if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
   qh_check_points();
 exitcode= qh_ERRnone;
 qh NOerrexit= True;  /* no more setjmp */
#ifdef qh_NOmem
  qh_freeqhull( True);
#else
  qh_freeqhull( False);
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
       totlong, curlong);
#endif
  return exitcode;
} /* main */


/****************************************************************************
* MeshLab                                                           o o     *
****************************************************************************/

//--- Include qhull, so it works from with in a C++ source file
//---
//--- In MVC one cannot just do:
//---
//---    extern "C"
//---    {
//---      #include "qhull_a.h"
//---    }
//---
//--- Because qhull_a.h includes math.h, which can not appear
//--- inside a extern "C" declaration.
//---
//--- Maybe that why Numerical recipes in C avoid this problem, by removing
//--- standard include headers from its header files and add them in the
//--- respective source files instead.
//---
//--- [K. Erleben]

/****************************************************************************
  History


****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <qhull.h>
#include <mem.h>
#include <qset.h>
#include <geom.h>
#include <merge.h>
#include <poly.h>
#include <io.h>
#include <stat.h>
#define false 0
#define true  1


facetT *compute_convex_hull(int dim, int numpoints,coordT * points);
facetT *compute_delaunay(int dim, int numpoints, coordT * points );
setT * compute_alpha_shapes(int dim, int numpoints, coordT * points, real alpha, int *nelements);
//real * qh_setradii(real * radii,int nradii);
//int visible_points(int dim, int numpoints, MeshModel &m, MeshModel &pm,MeshModel &pm2, vcg::Point3f viewpointP,float threshold,int convex_hullFP,int triangVP);

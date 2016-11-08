/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

/*
NOTE: how to acess atom properties
        if you have the index of the phase :
           for(i=0;i<itim->n[phase];i++){
                     atom_index = itim->gmx_index[phase][i];
                     resindex   = top->atoms.atom[atom_index].resind;
                     resname    = *top->atoms.resinfo[top->atoms.atom[atom_index].resind].name
                     atomname   = *(top->atoms.atomname[atom_index]);
                     pos_x      = itim->phase[phase][3*i];
           }
        if you have the index of the alpha-shape:
	   for(i=0;i<itim->nalphapoints;i++){
		    phase_index =itim->alpha_index[i];
                    atom_index  =itim->gmx_index[SUPPORT_PHASE][phase_index];
                    ...
           }

*/  

#if GMX_VERSION < 50000
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "string.h"
#include "smalloc.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "tpxio.h"
#include "gmx_ana.h"
#include <nbsearch.h>
#ifdef UNIX  
// used for handling signals
#include <stdio.h>
#include <unistd.h>
#include <sys/signal.h>
#endif //UNIX

#define gmx_ffopen ffopen
#define gmx_ffclose ffclose
#define wrap_gmx_rmpbc_init(a,b,c,d) gmx_rmpbc_init((a),(b),(c),(d))
#define wrap_read_next_x(a,b,c,d,e,f) read_next_x((a),(b),(c),(d),(e),(f))

#else

#define wrap_gmx_rmpbc_init(a,b,c,d) gmx_rmpbc_init((a),(b),(c))
#define wrap_read_next_x(a,b,c,d,e,f) read_next_x((a),(b),(c),(e),(f))
//#include "gmxpre.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlinehelpcontext.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "commandline/pargs.h"
#include "tpxio.h"
#include "trxio.h"
#include "xvgr.h"
#include "gmx_ana.h"
#include "gstat.h"
#include "macros.h"
#include "typedefs.h"
#include "viewit.h"
#include "vec.h"
#include "pbc.h"
#include "rmpbc.h"
#include "index.h"
#include "cstringutil.h"
#include "file.h"
#include "futil.h"
#include "smalloc.h"


#endif


#ifdef VIRIAL_EXTENSION
// this is defined in tpxio.h, in case...
#warning  COMPILING THE CODE USING THE VIRIAL_EXTENSION
#endif 


#ifdef TIME_PROFILE
#include <sys/time.h>
#endif

#define NORMAL_UNDEFINED -100
#define LAYER_OFFSET 10000

using namespace std;


int global_interrupt = 0;

typedef	enum { SUPPORT_PHASE=0, INNER_PHASE=1, OUTER_PHASE=2 } PHASE; // These value are not arbitrary and should not be changed.
 /* These value are not arbitrary and should not be changed: 
      OFF_NUMBER    -> for number density profile, used to normalize the order parameter density profile as well.
      OFF_ORDER[12] -> for order parameter density profile wrt macroscopic axis
      OFF_ORDER[34] -> for order parameter density profile wrt microscopic axis
  */
typedef	enum { OFF_NUMBER=4, OFF_ORDER1=7, OFF_ORDER1_2=10, OFF_ORDER2=13,  OFF_ORDER2_2=16, 
               OFF_ORDER3=19, OFF_ORDER3_2=22, OFF_ORDER4=25, OFF_ORDER4_2=28} HISTO_OFFSET; 
typedef enum {SURFACE_PLANE, SURFACE_SPHERE, SURFACE_CYLINDER, SURFACE_GENERIC } GEOMETRY;
typedef	enum { NONE, PATCH, FULL} PERIODIC ; 
typedef	enum { METHOD_ITIM, METHOD_A_SHAPE} METHOD; 

typedef struct {
	double  * rdata;
	double  size;
	double  minsize;
	int nbins;
	int N;
	int iterations;
	double bWidth;
	} Histogram;

typedef struct {
	struct kdtree * tree;
	int * flag;
	int n[2];
	real size[2];
 	int nelem;
} MESH;
/* TODO: put alphashape / surface variables in a separate structure, to be included here ?*/


typedef struct { 
	real box[3];
	real alpha;
	real * masses;
	real * charges;
	real skin;
	real range;
	real *pradii;
	int nphases;
	int RANDOM_PHASE;
	int ngmxphases;
	int normal;
	int info;
        int n_histo;
	int nadd_index;
	real ** phase;
 	int ** phase_index;
        int * inclusive_map;
        int dump_mol;
	int *n;
	real * alphapoints;
	int nalphapoints;
	int current_layer;
	int maxlayers;
	int * alpha_index;
	real * radii;
        int * mask;
	int ** mask_add;
	int * gmx_alpha_id;
	int * gmx_alpha_phase_id;
	int ** gmx_index;
	int *indexm; // for the molecular version of the support group.
	int *backindex; // for the molecular version of the support group.
	int com_opt[64];
        int bCom;
        int bVirial;
        int bMol;
        int bInclusive;
        int bMCnormalization;
        int bOrder;
	int side;
	real target_mesh_size;
	MESH mesh;
	GEOMETRY geometry;
	PERIODIC *periodic;
	METHOD method;
	char method_name[2][128];
	void (* dump_surface_points)(t_topology*,FILE*);		
	void (* dump_slabs)(t_topology*);		
	void (* dump_surface_molecules)(t_topology*,FILE*,atom_id **);		
	void (* dump_phase_points)(PHASE,t_topology*,FILE*);		
	Histogram * histograms;

} ITIM;

ITIM * global_itim;

#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

#define my_min(a,b) ( (a) > (b) ? (b) : (a))
double *CLUSTERCUT;     /* group dependent cut-off for cluster analysis */ // gh 
double *CLUSTERCUT2;     /* group dependent cut-off for cluster analysis */ // gh 

// TODO: make it an option? change the algorithm ? 
#define MAX_KDTREE_CHECK 60

struct kdtree;
struct kdres;

//typedef float real;
/* create a kd-tree for "k"-dimensional data */
struct kdtree *kd_create(int k);

/* free the struct kdtree */
void kd_free(struct kdtree *tree);

/* remove all the elements from the tree */
void kd_clear(struct kdtree *tree);

/* if called with non-null 2nd argument, the function provided
 * will be called on data pointers (see kd_insert) when nodes
 * are to be removed from the tree.
 */
void kd_data_destructor(struct kdtree *tree, void (*destr)(void*));

/* insert a node, specifying its position, and optional data */
int kd_insert(struct kdtree *tree, const real *pos, real *data);
int kd_insertf(struct kdtree *tree, const float *pos, real *data);
int kd_insert3(struct kdtree *tree, real x, real y, real z, real *data);
int kd_insert3p(struct kdtree *tree, real x, real y, real z, real *data,real cut);
int kd_insert3f(struct kdtree *tree, float x, float y, float z, real *data);

/* Find one of the nearest nodes from the specified point.
 *
 * This function returns a pointer to a result set with at most one element.
 */
struct kdres *kd_nearest(struct kdtree *tree, const real *pos);
struct kdres *kd_nearestf(struct kdtree *tree, const float *pos);
struct kdres *kd_nearest3(struct kdtree *tree, real x, real y, real z);
struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z);

/* Find any nearest nodes from the specified point within a range.
 *
 * This function returns a pointer to a result set, which can be manipulated
 * by the kd_res_* functions.
 * The returned pointer can be null as an indication of an error. Otherwise
 * a valid result set is always returned which may contain 0 or more elements.
 * The result set must be deallocated with kd_res_free, after use.
 */
struct kdres *kd_nearest_range(struct kdtree *tree, const real *pos, real range);
struct kdres *kd_nearest_rangef(struct kdtree *tree, const float *pos, float range);
struct kdres *kd_nearest_range3(struct kdtree *tree, real x, real y, real z, real range);
struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range);

/* frees a result set returned by kd_nearest_range() */
void kd_res_free(struct kdres *set);

/* returns the size of the result set (in elements) */
int kd_res_size(struct kdres *set);

/* rewinds the result set iterator */
void kd_res_rewind(struct kdres *set);

/* returns non-zero if the set iterator reached the end after the last element */
int kd_res_end(struct kdres *set);

/* advances the result set iterator, returns non-zero on success, zero if
 * there are no more elements in the result set.
 */
int kd_res_next(struct kdres *set);

/* returns the data pointer (can be null) of the current result set item
 * and optionally sets its position to the pointers(s) if not null.
 */
void *kd_res_item(struct kdres *set, real *pos);
void *kd_res_itemf(struct kdres *set, float *pos);
void *kd_res_item3(struct kdres *set, real *x, real *y, real *z);
void *kd_res_item3f(struct kdres *set, float *x, float *y, float *z);

/* equivalent to kd_res_item(set, 0) */
void *kd_res_item_data(struct kdres *set);

real compute_osculating_sphere_radius(real p[3], real q[3], real r[3] ,real s[3], real wp, real wq, real wr, real ws);

real interpolate_distance3D(real *A,real *B,real *C,real *I);

#ifdef __cplusplus
}
#endif

#endif	/* _KDTREE_H_ */
#define NO_ALLOCA
#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#ifdef USE_LIST_NODE_ALLOCATOR

#ifndef NO_PTHREADS
#include <pthread.h>
#else

#ifndef I_WANT_THREAD_BUGS
#error "You are compiling with the fast list node allocator, with pthreads disabled! This WILL break if used from multiple threads."
#endif	/* I want thread bugs */

#endif	/* pthread support */
#endif	/* use list node allocator */

struct kdhyperrect {
	int dim;
	real *min, *max;              /* minimum/maximum coords */
};

struct kdnode {
	real *pos;
	int dir;
	real *data;

	struct kdnode *left, *right;	/* negative/positive side */
};

struct res_node {
	struct kdnode *item;
	real dist_sq;
	struct res_node *next;
};

struct kdtree {
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
	void (*destr)(void*);
};

struct kdres {
	struct kdtree *tree;
	struct res_node *rlist, *riter;
	int size;
};

#define SQ(x)			((x) * (x))


static void clear_rec(struct kdnode *node, void (*destr)(void*));
static int insert_rec(struct kdnode **node, const real *pos, real *data, int dir, int dim);
static int rlist_insert(struct res_node *list, struct kdnode *item, real dist_sq);
static void clear_results(struct kdres *set);

static struct kdhyperrect* hyperrect_create(int dim, const real *min, const real *max);
static void hyperrect_free(struct kdhyperrect *rect);
static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect);
static void hyperrect_extend(struct kdhyperrect *rect, const real *pos);
static real hyperrect_dist_sq(struct kdhyperrect *rect, const real *pos);

#ifdef USE_LIST_NODE_ALLOCATOR
static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node*);
#else
#define alloc_resnode()		(struct res_node*)malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif



struct kdtree *kd_create(int k)
{
	struct kdtree *tree;

	if(!(tree = (struct kdtree*)malloc(sizeof *tree))) {
		return 0;
	}

	tree->dim = k;
	tree->root = 0;
	tree->destr = free;
	tree->rect = 0;

	return tree;
}

void kd_free(struct kdtree *tree)
{
	if(tree) {
		kd_clear(tree);
		free(tree);
	}
}

static void clear_rec(struct kdnode *node, void (*destr)(void*))
{
	if(!node) return;

	clear_rec(node->left, destr);
	clear_rec(node->right, destr);
	
//	if(destr) {
//		destr(node->data);
//	}
	free(node->pos);
	free(node);
}

void kd_clear(struct kdtree *tree)
{
	clear_rec(tree->root, tree->destr);
	tree->root = 0;

	if (tree->rect) {
		hyperrect_free(tree->rect);
		tree->rect = 0;
	}
}

void kd_data_destructor(struct kdtree *tree, void (*destr)(void*))
{
	tree->destr = destr;
}


static int insert_rec(struct kdnode **nptr, const real *pos, real *data, int dir, int dim)
{
	int new_dir;
	struct kdnode *node;

	if(!*nptr) {
		if(!(node = (struct kdnode*)malloc(sizeof *node))) {
			return -1;
		}
		if(!(node->pos = (real*)malloc(dim * sizeof *node->pos))) {
			free(node);
			return -1;
		}
//		if(!(node->data = (real *) malloc(sizeof *node->data))) {
//			free(node);
//			return -1;
//		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
//		memcpy(node->data, data,  sizeof *node->data);
		node->data = data;
		node->dir = dir;
		node->left = node->right = 0;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if(pos[node->dir] < node->pos[node->dir]) {
		return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
	}
	return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

int kd_insert(struct kdtree *tree, const real *pos, real *data)
{
	if (insert_rec(&tree->root, pos, data, 0, tree->dim)) {
		return -1;
	}

	if (tree->rect == 0) {
		tree->rect = hyperrect_create(tree->dim, pos, pos);
	} else {
		hyperrect_extend(tree->rect, pos);
	}

	return 0;
}

int kd_insertf(struct kdtree *tree, const float *pos, real *data)
{
	static real sbuf[16];
	real *bptr, *buf = 0;
	int res, dim = tree->dim;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (real*)malloc(dim * sizeof *bptr))) {
				return -1;
			}
	} else {
		bptr = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_insert(tree, buf, data);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
		free(buf);
	return res;
}

int kd_insert3(struct kdtree *tree, real x, real y, real z, real *data)
{
	real buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

int kd_insert3p(struct kdtree *tree, real x, real y, real z, real *data, real cut)
{
	ITIM* itim = global_itim;
	int i,j,k,ret,count=0;
	real buf[3],buf2[3];
	buf[0]=x;
	buf[1]=y;
	buf[2]=z;
	ret = kd_insert(tree, buf, data);
	if(cut>0){
	     for(i=-1;i<2;i++){
	       if( (i==-1 && buf[0]+cut > itim->box[0]/2.) || (i==1 && buf[0]-cut < -itim->box[0]/2.) || i==0 ){
	         buf2[0]=buf[0]+i*itim->box[0];	
	         for(j=-1;j<2;j++){
	           if( (j==-1 && buf[1]+cut > itim->box[1]/2.) || (j==1 && buf[1]-cut < -itim->box[1]/2.) || j==0 ){
	             buf2[1]=buf[1]+j*itim->box[1];	
	             for(k=-1;k<2;k++){
	                if( (k==-1 && buf[2]+cut > itim->box[2]/2.) || (k==1 && buf[2]-cut < -itim->box[2]/2.) || k==0 ){
	                  buf2[2]=buf[2]+k*itim->box[2];	
	                  if(k==0 && j==0 && i==0) continue;
   	     		  ret = kd_insert(tree, buf2, data);
			  count++;
	                }
	             }
                   }
                 }
	       }
             }
        }
	return ret;
}

int kd_insert3f(struct kdtree *tree, float x, float y, float z, real *data)
{
	real buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

static int find_nearest(struct kdnode *node, const real *pos, real range, struct res_node *list, int ordered, int dim)
{
	real dist_sq, dx;
	int i, ret, added_res = 0;
	ITIM * itim = global_itim ;
	if(!node) return 0;

	dist_sq = 0;
        if(itim!=NULL) { 
	  for(i=0; i<dim; i++) {
	  	dist_sq += SQ(node->pos[i] - pos[i]);
	  }
        } else { 
	  for(i=0; i<dim; i++) {
		real  dist = node->pos[i] - pos[i];
		while(dist > itim->box[i]/2.)  dist-=itim->box[i];
		while(dist < -itim->box[i]/2.) dist+=itim->box[i];
	  	dist_sq += SQ(dist);
	  }
        }
	if(dist_sq <= SQ(range)) {
		if(rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1) {
			return -1;
		}
		added_res = 1;
	}

	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest(dx <= 0.0 ? node->left : node->right, pos, range, list, ordered, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest(dx <= 0.0 ? node->right : node->left, pos, range, list, ordered, dim);
	}
	if(ret == -1) {
		return -1;
	}
	added_res += ret;

	return added_res;
}

static void kd_nearest_i(struct kdnode *node, const real *pos, struct kdnode **result, real *result_dist_sq, struct kdhyperrect* rect)
{
	int dir = node->dir;
	int i, side;
	real dummy, dist_sq;
	struct kdnode *nearer_subtree, *farther_subtree;
	real *nearer_hyperrect_coord, *farther_hyperrect_coord;

	/* Decide whether to go left or right in the tree */
	dummy = pos[dir] - node->pos[dir];
	if (dummy <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = rect->max + dir;
		farther_hyperrect_coord = rect->min + dir;
		side = 0;
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = rect->min + dir;
		farther_hyperrect_coord = rect->max + dir;
		side = 1;
	}

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		kd_nearest_i(nearer_subtree, pos, result, result_dist_sq, rect);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

	/* Check the distance of the point at the current node, compare it
	 * with our best so far */
	dist_sq = 0;
	for(i=0; i < rect->dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if (dist_sq < *result_dist_sq) {
		*result = node;
		*result_dist_sq = dist_sq;
	}

	if (farther_subtree) {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our
		 * minimum distance in result_dist_sq. */
		if (hyperrect_dist_sq(rect, pos) < *result_dist_sq) {
			/* Recurse down into farther subtree */
			kd_nearest_i(farther_subtree, pos, result, result_dist_sq, rect);
		}
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
	}
}

struct kdres *kd_nearest(struct kdtree *kd, const real *pos)
{
	struct kdhyperrect *rect;
	struct kdnode *result;
	struct kdres *rset;
	real dist_sq;
	int i;

	if (!kd) return 0;
	if (!kd->rect) return 0;

	/* Allocate result set */
	if(!(rset = (struct kdres *)malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	/* Duplicate the bounding hyperrectangle, we will work on the copy */
	if (!(rect = hyperrect_duplicate(kd->rect))) {
		kd_res_free(rset);
		return 0;
	}

	/* Our first guesstimate is the root node */
	result = kd->root;
	dist_sq = 0;
	for (i = 0; i < kd->dim; i++)
		dist_sq += SQ(result->pos[i] - pos[i]);

	/* Search for the nearest neighbour recursively */
	kd_nearest_i(kd->root, pos, &result, &dist_sq, rect);

	/* Free the copy of the hyperrect */
	hyperrect_free(rect);

	/* Store the result */
	if (result) {
		if (rlist_insert(rset->rlist, result, -1.0) == -1) {
			kd_res_free(rset);
			return 0;
		}
		rset->size = 1;
		kd_res_rewind(rset);
		return rset;
	} else {
		kd_res_free(rset);
		return 0;
	}
}

struct kdres *kd_nearestf(struct kdtree *tree, const float *pos)
{
	static real sbuf[16];
	real *bptr, *buf = 0;
	int dim = tree->dim;
	struct kdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (real *)malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_nearest(tree, buf);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
		free(buf);
	return res;
}

struct kdres *kd_nearest3(struct kdtree *tree, real x, real y, real z)
{
	real pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return kd_nearest(tree, pos);
}

struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z)
{
	real pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return kd_nearest(tree, pos);
}

struct kdres *kd_nearest_range(struct kdtree *kd, const real *pos, real range)
{
	int ret;
	struct kdres *rset;
        
	if(!(rset = (struct kdres *) malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	if((ret = find_nearest(kd->root, pos, range, rset->rlist, 1, kd->dim)) == -1) {
		kd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	kd_res_rewind(rset);
	return rset;
}

struct kdres *kd_nearest_rangef(struct kdtree *kd, const float *pos, float range)
{
	static real sbuf[16];
	real *bptr, *buf = 0;
	int dim = kd->dim;
	struct kdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (real *)malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_nearest_range(kd, buf, range);
#ifndef NO_ALLOCA
	if(kd->dim > 256)
#else
	if(kd->dim > 16)
#endif
		free(buf);
	return res;
}

struct kdres *kd_nearest_range3(struct kdtree *tree, real x, real y, real z, real range)
{
	real buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}

struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range)
{
	real buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}

void kd_res_free(struct kdres *rset)
{
	clear_results(rset);
	free_resnode(rset->rlist);
	free(rset);
}

int kd_res_size(struct kdres *set)
{
	return (set->size);
}

void kd_res_rewind(struct kdres *rset)
{
	rset->riter = rset->rlist->next;
}

int kd_res_end(struct kdres *rset)
{
	return rset->riter == 0;
}

int kd_res_next(struct kdres *rset)
{
	rset->riter = rset->riter->next;
	return rset->riter != 0;
}

void *kd_res_item(struct kdres *rset, real *pos)
{
	if(rset->riter) {
		if(pos) {
			memcpy(pos, rset->riter->item->pos, rset->tree->dim * sizeof *pos);
		}
		return rset->riter->item->data;
	}
	return 0;
}

void *kd_res_itemf(struct kdres *rset, float *pos)
{
	if(rset->riter) {
		if(pos) {
			int i;
			for(i=0; i<rset->tree->dim; i++) {
				pos[i] = rset->riter->item->pos[i];
			}
		}
		return rset->riter->item->data;
	}
	return 0;
}

void *kd_res_item3(struct kdres *rset, real *x, real *y, real *z)
{
	if(rset->riter) {
		if(*x) *x = rset->riter->item->pos[0];
		if(*y) *y = rset->riter->item->pos[1];
		if(*z) *z = rset->riter->item->pos[2];
	}
	return 0;
}

void *kd_res_item3f(struct kdres *rset, float *x, float *y, float *z)
{
	if(rset->riter) {
		if(*x) *x = rset->riter->item->pos[0];
		if(*y) *y = rset->riter->item->pos[1];
		if(*z) *z = rset->riter->item->pos[2];
	}
	return 0;
}

void *kd_res_item_data(struct kdres *set)
{
	return kd_res_item(set, 0);
}

/* ---- hyperrectangle helpers ---- */
static struct kdhyperrect* hyperrect_create(int dim, const real *min, const real *max)
{
	size_t size = dim * sizeof(real);
	struct kdhyperrect* rect = 0;

	if (!(rect =(struct kdhyperrect* )malloc(sizeof(struct kdhyperrect)))) {
		return 0;
	}

	rect->dim = dim;
	if (!(rect->min = (real *)malloc(size*sizeof(real*)))) {
		free(rect);
		return 0;
	}
	if (!(rect->max = (real *)malloc(size*sizeof(real*)))) {
		free(rect->min);
		free(rect);
		return 0;
	}
	memcpy(rect->min, min, size);
	memcpy(rect->max, max, size);

	return rect;
}

static void hyperrect_free(struct kdhyperrect *rect)
{
	free(rect->min);
	free(rect->max);
	free(rect);
}

static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect)
{
	return hyperrect_create(rect->dim, rect->min, rect->max);
}

static void hyperrect_extend(struct kdhyperrect *rect, const real *pos)
{
	int i;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			rect->min[i] = pos[i];
		}
		if (pos[i] > rect->max[i]) {
			rect->max[i] = pos[i];
		}
	}
}

static real hyperrect_dist_sq(struct kdhyperrect *rect, const real *pos)
{
	int i;
	real result = 0;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			result += SQ(rect->min[i] - pos[i]);
		} else if (pos[i] > rect->max[i]) {
			result += SQ(rect->max[i] - pos[i]);
		}
	}

	return result;
}

/* ---- static helpers ---- */

#ifdef USE_LIST_NODE_ALLOCATOR
/* special list node allocators. */
static struct res_node *free_nodes;

#ifndef NO_PTHREADS
static pthread_mutex_t alloc_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

static struct res_node *alloc_resnode(void)
{
	struct res_node *node;

#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	if(!free_nodes) {
		node = (struct res_node *)malloc(sizeof *node);
	} else {
		node = free_nodes;
		free_nodes = free_nodes->next;
		node->next = 0;
	}

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif

	return node;
}

static void free_resnode(struct res_node *node)
{
#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	node->next = free_nodes;
	free_nodes = node;

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif
}
#endif	/* list node allocator or not */


/* inserts the item. if dist_sq is >= 0, then do an ordered insert */
static int rlist_insert(struct res_node *list, struct kdnode *item, real dist_sq)
{
	struct res_node *rnode;

	if(!(rnode = alloc_resnode())) {
		return -1;
	}
	rnode->item = item;
	rnode->dist_sq = dist_sq;

	if(dist_sq >= 0.0) {
		while(list->next && list->next->dist_sq < dist_sq) {
			list = list->next;
		}
	}
	rnode->next = list->next;
	list->next = rnode;
	return 0;
}

static void clear_results(struct kdres *rset)
{
	struct res_node *tmp, *node = rset->rlist->next;

	while(node) {
		tmp = node;
		node = node->next;
		free_resnode(tmp);
	}

	rset->rlist->next = 0;
}



/* NOTE: alpha shape code adapted from meshlab */
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#define HISTO_N 100

#define SQR(x) ((x)*(x))



#define determ2d( a1,a2,b1,b2 ) (( a1 )*( b2 ) - ( a2 )*( b1 ))

#define determ3d( a1,a2,a3,b1,b2,b3,c1,c2,c3 ) ( ( a1 )*determ2d( b2,b3,c2,c3 ) \
                - ( b1 )*determ2d( a2,a3,c2,c3 ) + ( c1 )*determ2d( a2,a3,b2,b3 ) )

/* TODO: These should disappear when including functions into g_density */
//typedef real  matrix[3][3];
//typedef int atom_id;
/**/

// to be removed when linked to gromacs

//typedef real rvec[3];
//
typedef enum { DIR_NEGATIVE=-1, DIR_POSITIVE=1} Direction;
struct kdtree *Surface=NULL;
struct kdtree *Bulk=NULL;
typedef struct { 
	real box[3];
	real ref[3];
} qsort_data;

real * global_real_pointer;
real * global_masses=NULL;
real * global_charges=NULL;

typedef struct {
   int bCom ; 
   int bMol ;
   int bVirial ;
   int bMCnormalization ; 
   char com_opt_file[1024];
   gmx_bool bCenter ;
   gmx_bool bInfo ;
   gmx_bool bIntrinsic ;
   gmx_bool bInclusive ;
   gmx_bool bCluster ;
   gmx_bool bDump ;
   gmx_bool bDumpPhases ;
   gmx_bool bOrder ;
   gmx_bool bHelp ;
   gmx_bool bSymmetrize;
 
} Flags ; 

void InitFlags(Flags *f){
  f->bCom=f->bMol=f->bVirial=FALSE;
  f->bMCnormalization=0; 
  f->bCenter=f->bInfo=f->bInclusive=f->bCluster=FALSE;
  f->bDump=f->bDumpPhases=f->bOrder=f->bHelp=f->bSymmetrize=FALSE;
  f->bIntrinsic=TRUE;
}





#ifdef INTEGER_HISTO
unsigned long long int * int_histo;
#endif
ITIM * init_intrinsic_surface(Flags myflags, int normal, real alpha, real mesh,  matrix box, int ngrps,int *nbins, int maxlayers,real * radii, int ** index, int *gnx,int *com_opt, int dump_mol, const char ** geometry, int ngrps_add, t_topology * top);


void  compute_intrinsic_profile(matrix box, atom_id **index, t_topology * top, char dens_opt, t_trxframe * fr);
void  compute_intrinsic_order(matrix box, atom_id **index, t_topology * top);
void  finalize_intrinsic_profile(real *** density, int * nslices, real * slWidth);
void  free_profile(real *** density);


int add_histo(Histogram * histo, int N, double pos, double value){
	int bin ;
	bin = (int) (pos / histo->bWidth);
        if(bin<histo->nbins && bin >=0){ 
               histo->rdata[N*histo->nbins+bin] += value;
#ifdef INTEGER_HISTO
               int_histo[N*histo->nbins+bin]++;
#endif
               return 1 ;
        } 
        return 0;
}
void dump_histo(Histogram * histo,int N, FILE*file){
	int i;
	for(i=0;i<histo->nbins;i++) fprintf(file,"%f %f\n",i*histo->bWidth,histo->rdata[N*histo->nbins+i]*1.0/histo->iterations);
}

typedef enum {LAYER_DISTRIBUTION=0,INTRINSIC_DENSITY=1,INTRINSIC_ORDER=2,COMPUTE_SIZE=3} HISTO_TYPE; 
typedef enum {ATOMIC=0,MOLECULAR=1} HISTO_MOLTYPE; 

int GET_HISTO_INDEX(HISTO_TYPE TYPE,int phase, int layer, int molecular, int line){
	int index = 0, i,moltype,n;
	ITIM * itim = global_itim;
	for (moltype=ATOMIC;moltype<=MOLECULAR;moltype++){
      		if (itim->bMol==FALSE && moltype==MOLECULAR) continue;
        	if(layer!=0 && !(phase==SUPPORT_PHASE || phase >= itim->ngmxphases)) { 
        	    exit(printf("Internal error: layer intrinsic profile only possible for the support phase %d %d\n",phase,itim->ngmxphases));
        	}
        	for(i=1;i<itim->maxlayers+1;i++){ // two (layer + intrinsic) histograms x number of layers for the support phase 
        	    	if(moltype == molecular && TYPE==LAYER_DISTRIBUTION && layer==i  && phase==SUPPORT_PHASE) return index;
        	    	index++;
        	    	if(moltype == molecular && TYPE==INTRINSIC_DENSITY  && layer==i  && phase==SUPPORT_PHASE) return index;
        	    	index++;
        	}
        }
	for(int add=0;add<itim->nadd_index;add++){
	  for (moltype=ATOMIC;moltype<=MOLECULAR;moltype++){
      		if (itim->bMol==FALSE && moltype==MOLECULAR) continue;
          	if(layer!=0 && !(phase==SUPPORT_PHASE || phase >= itim->ngmxphases)) { 
          	    exit(printf("Internal error: layer intrinsic profile only possible for the support phase\n"));
          	}
          	for(i=1;i<itim->maxlayers+1;i++){ // two (layer + intrinsic) histograms x number of layers for the support phase 
          	    	if(moltype == molecular && TYPE==LAYER_DISTRIBUTION && layer==i  && phase==itim->ngmxphases+add) return index;
          	    	index++;
          	    	if(moltype == molecular && TYPE==INTRINSIC_DENSITY  && layer==i  && phase==itim->ngmxphases+add) return index;
          	    	index++;
          	}
          }
        }
        for(n=SUPPORT_PHASE+1;n<itim->nphases;n++){ // two (nonintrinsic + intrinsic) histogram for each of the other phases
            	if(molecular == ATOMIC && phase==n && TYPE==LAYER_DISTRIBUTION) return index;
            	index++;
            	if(molecular == ATOMIC && phase==n && TYPE==INTRINSIC_DENSITY) return index;
            	index++;
        }
        if (TYPE==INTRINSIC_ORDER) {
                exit(printf("Internal error: order profiles to be checked\n"));
        } 
       if(TYPE==COMPUTE_SIZE){
       	return index;   	
       } else { 
	   return -1000;
	// exit(printf("Internal error: GET_HISTO_INDEX() called with wrong arguments at line %d: TYPE=%d, phase=%d, molecular=%d, layer=%d\n",line,TYPE,phase,molecular,layer));
       }
}

int GET_HISTO_SIZE(void) {
 return GET_HISTO_INDEX(COMPUTE_SIZE,0, 0, 0,__LINE__);
}

void plot_intrinsic_density(Histogram * histo,char **grpmname, const char * fn, char dens_opt){
   ITIM*itim=global_itim;
   int nphases = itim->nphases-1;
   FILE *cid = fopen(fn,"w");
   static const char *modif[]={"atomic","molecular"} ;
   int column=1,moltype,i,j,n;
   real factor=1;
   switch(dens_opt){
	case 'm': factor = AMU/(NANO*NANO*NANO); break;
   }
   fprintf(cid,"#column %d : position \n",column); column++;		
   for (moltype=0;moltype<=MOLECULAR-ATOMIC;moltype++){
      if (itim->bMol==FALSE && moltype==MOLECULAR) continue;
      for(i=1;i<itim->maxlayers+1;i++){
	fprintf(cid,"#column %d : support group (%s) nonintr. %s dens. layer %d\n",column,grpmname[0],modif[moltype],i); column++;		
	fprintf(cid,"#column %d : support group (%s) intrins. %s dens. layer %d\n",column,grpmname[0],modif[moltype],i); column++;		
      }
   }
   for(int add=0;add<itim->nadd_index;add++){
     for (moltype=0;moltype<=MOLECULAR-ATOMIC;moltype++){
        if (itim->bMol==FALSE && moltype==MOLECULAR) continue;
        for(i=1;i<itim->maxlayers+1;i++){
          fprintf(cid,"#column %d : additional group (%s)  nonintr. dens. %s layer %d\n",column,grpmname[itim->ngmxphases+add],modif[moltype],i); column++;		
          fprintf(cid,"#column %d : additional group (%s)  intrins. dens. %s layer %d\n",column,grpmname[itim->ngmxphases+add],modif[moltype],i); column++;		
        }
     }
   }

   for(n=SUPPORT_PHASE+1;n<nphases;n++){
     fprintf(cid,"#column %d : group %d (%s)  nonintr. dens. %s\n",column,n-SUPPORT_PHASE,grpmname[n],modif[ATOMIC]); column++;
     fprintf(cid,"#column %d : group %d (%s)  intrins. dens. %s\n",column,n-SUPPORT_PHASE,grpmname[n],modif[ATOMIC]); column++;
     fflush(cid);
   }
   for(j=0;j<histo->nbins;j++) { 
     fprintf(cid,"%f ",(j-histo->nbins/2)*histo->bWidth);
     for (i=0;i<itim->n_histo -2 ;i++){ // n_histo -2 to skip the random phase.
        fprintf(cid," %f ",histo->rdata[i*histo->nbins+j] * factor);
     }
     fprintf(cid,"\n");
   }
}

#ifdef INTEGER_HISTO
void dump_int_histo(Histogram * histo,int N, FILE*file){
	int i;
        fprintf(file,"# iterations = %d   normalization to number density= %f\n",histo->iterations,4*(global_itim->box[0] * global_itim->box[1])*histo->bWidth );
	for(i=0;i<histo->nbins;i++) fprintf(file,"%f %llu\n",i*histo->bWidth,int_histo[N*histo->nbins+i]);
}
#endif
void clear_histo(Histogram * histo,int N){
	int i;
	for(i=0;i<histo->nbins;i++) histo->rdata[N*histo->nbins+i]=0;
}
void eprintv(real *v, char * s){
	fprintf(stdout,"%f %f %f",v[0],v[1],v[2]);
	fprintf(stdout,"%s",s);
}
void printv(real *v, char * s){
	printf("%f %f %f",v[0],v[1],v[2]);
	printf("%s",s);
}
real remove_pbc_1d(real a, real box){
	while(a>box/2.) a-=box; 
	while(a<-box/2.) a+=box; 
	return a;
}

real triangle_area(real *p0, real *p1, real *p2,real *box)
{
	int i;
	real d1[2],d2[2];
        for(i=0;i<2;i++){
           d1[i]=p1[i]-p0[i];
           d2[i]=p2[i]-p0[i];
	   while (d1[i]>box[i]/2.) d1[i]-=box[i];
	   while (d1[i]<-box[i]/2.) d1[i]+=box[i];
	   while (d2[i]>box[i]/2.) d2[i]-=box[i];
	   while (d2[i]<-box[i]/2.) d2[i]+=box[i];
	}
	return 0.5*( d1[0]*d2[1] - d2[0]*d1[1] );
}





real interpolate_distance(real *A,real *B,real *C,real *I){
	int i;
	real ba[2],cb[2],bi[2],ci[2],ai[2],ac[2];
        for(i=0;i<2;i++){
		ba[i]=remove_pbc_1d( B[i]-A[i]  ,global_itim->box[i]);
		cb[i]=remove_pbc_1d( C[i]-B[i]  ,global_itim->box[i]);
		ac[i]=remove_pbc_1d( A[i]-C[i]  ,global_itim->box[i]);
		bi[i]=remove_pbc_1d( B[i]-I[i]  ,global_itim->box[i]);
		ci[i]=remove_pbc_1d( C[i]-I[i]  ,global_itim->box[i]);
		ai[i]=remove_pbc_1d( A[i]-I[i]  ,global_itim->box[i]);
	}
	real  F = ba[0] * cb[1] - ba[1] * cb[0] ;
	return A[2]*(  bi[0]*cb[1] - bi[1]*cb[0] )  / F +
	       B[2]*(  ci[0]*ac[1] - ci[1]*ac[0] )  / F +
	       C[2]*(  ai[0]*ba[1] - ai[1]*ba[0] )  / F ;
}

void arrange_bulk_points (int Bulk_nelem, real * Bulk_points, GEOMETRY geometry, int dimension ){ 
	int i;
	Bulk  = kd_create( dimension );
	for(i=0;i<Bulk_nelem;i++){
		switch(dimension) { 
			case 2: 
			        kd_insert(Bulk, &Bulk_points[3*i], &Bulk_points[3*i+2]); 
			    break;
			case 3:
			       kd_insert(Bulk,&Bulk_points[3*i],0 ); 
			    break;
		}

	}
}

void collect_statistics_for_layers(ITIM * itim, t_trxframe * fr){
	int i,layer,atom_index,phase_index;
	real * f,*vir,*p,*J;
	real Jtot[3];
	int *n;
	FILE*fp;
	n   = (int*)malloc(itim->maxlayers * sizeof(int));
	p   = (real*)malloc(itim->maxlayers * sizeof(real));
	f   = (real*)malloc(itim->maxlayers * 3 * sizeof(real));
	vir = (real*)malloc(itim->maxlayers * 3 * sizeof(real));
	J   = (real*)malloc((itim->maxlayers+1) * 3 * sizeof(real));
	for(i = 0 ; i < itim->maxlayers ; i++ ){
		n[i]=0;
		p[i]=0;
		f[3*i]=f[3*i+1]=f[3*i+2]=0.0;
		vir[3*i]=vir[3*i+1]=vir[3*i+2]=0.0;
		J[3*i]=J[3*i+1]=J[3*i+2]=0.0;
		J[3*(i+1)]=J[3*(i+1)+1]=J[3*(i+1)+2]=0.0;
	}
	Jtot[0]=Jtot[1]=Jtot[2]=0.0;
 	for(i=0;i<itim->n[SUPPORT_PHASE];i++){
	       atom_index = itim->gmx_index[SUPPORT_PHASE][i];
               double m = itim->masses[atom_index];
#ifdef VIRIAL_EXTENSION
	       Jtot[0] += fr->v[atom_index][0] * itim->charges[atom_index];
	       Jtot[1] += fr->v[atom_index][1] * itim->charges[atom_index];
	       Jtot[2] += fr->v[atom_index][2] * itim->charges[atom_index];
#endif
	       if(itim->mask[i]>0 && itim->mask[i]< itim->maxlayers+1){
		 layer = itim->mask[i]-1;
		 n[layer]+=1;
		 real sign=(itim->phase[SUPPORT_PHASE][3*i+itim->normal]>0?1:-1);
                 p[layer]+=sign*itim->phase[SUPPORT_PHASE][3*i+itim->normal]*itim->charges[atom_index];
#ifdef VIRIAL_EXTENSION
	         J[3*layer+0] += fr->v[atom_index][0] * itim->charges[atom_index];
	         J[3*layer+1] += fr->v[atom_index][1] * itim->charges[atom_index];
	         J[3*layer+2] += fr->v[atom_index][2] * itim->charges[atom_index];
                 f[3*layer]   += fr->f[atom_index][0];
                 f[3*layer+1] += fr->f[atom_index][1];
                 f[3*layer+2] += fr->f[atom_index][2];
                 vir[3*layer]   +=  fr->vir[atom_index][0];
                 vir[3*layer+1] +=  fr->vir[atom_index][1];
                 vir[3*layer+2] +=  fr->vir[atom_index][2];
#endif
               }
        }
	fp = fopen("stats.dat","a");	

	real surface = itim->box[0] * itim->box[1];
	real volume =  surface * itim->box[2];
	for(i = 0 ; i < itim->maxlayers ; i++ ){
		fprintf(fp, "%f ",n[i]/surface);
		fprintf(fp, "%f ",p[i]/n[i]);
		fprintf(fp,"%f %f %f ", J[3*i],J[3*i+1],J[3*i+2]);
#ifdef VIRIAL_EXTENSION
		fprintf(fp, "%f %f %f ",vir[3*i]/volume,vir[3*i+1]/volume,vir[3*i+2]/volume);
#endif
        }
	fprintf(fp,"%f %f %f ", Jtot[0],Jtot[1],Jtot[2]);
	fprintf(fp, "\n");
	fclose(fp);
	free(n);
	free(p);
	free(f);
	free(J);
	free(vir);
}

void init_itim_grid(ITIM * itim){
	int i;
	real pos[2];
	for(i=0;i<2;i++){
   	  itim->mesh.n[i] = ceil(itim->box[i]/itim->target_mesh_size);
   	  itim->mesh.size[i] = itim->box[i]/itim->mesh.n[i];
	}
	itim->mesh.nelem = itim->mesh.n[0]*itim->mesh.n[1];
        itim->mesh.flag = (int *) realloc (itim->mesh.flag,itim->mesh.nelem * sizeof(int));
	if( itim->mesh.tree != NULL) {
			kd_free( itim->mesh.tree);
			itim->mesh.tree=NULL;
	}
	itim->mesh.tree = kd_create(2);
	pos[0] = - itim->box[0]/2 ; 
	pos[1] = - itim->box[1]/2 ;
	for(i=0;i<itim->mesh.nelem;i++){
             	itim->mesh.flag[i]=DIR_POSITIVE;   
		kd_insert(itim->mesh.tree, pos, (real*)&(itim->mesh.flag[i])); 
		pos[0] += itim->mesh.size[0]; 
		if(pos[0] >= itim->box[0]/2){
			pos[1]+= itim->mesh.size[1]; 
			pos[0] =  - itim->box[0]/2 ;
		}
	}
}
void dump_itim_grid(ITIM * itim){
	struct kdres *presults;
	real pos[2]={0.0,0.0};
	real res[2];
	int *p;
	presults = kd_nearest_range( itim->mesh.tree, pos, 10*itim->box[0]); 
        while( !kd_res_end( presults ) ) {
          p =  (int*)kd_res_item(presults, res);
		printf("mesh %f %f %d\n",res[0],res[1],*p);
          /* go to the next entry */
          kd_res_next( presults );
        }
	kd_res_free( presults );
}

int is_in_correct_phase(int index){
	ITIM * itim=global_itim;
	int test=0;
	if(itim->bInclusive) {
			if( index>=itim->inclusive_map[2*INNER_PHASE] && index<=itim->inclusive_map[2*INNER_PHASE+1]) {
				 if( itim->mask[index]==0 || itim->mask[index] >= LAYER_OFFSET) {  test=1; }
			} else { 
				 if( itim->mask[index]==-1) {  test=1; }
			}
		
	} else { test = (itim->mask[index]==0 || itim->mask[index] >= LAYER_OFFSET) ; } 
	return test;
}

int check_itim_testlines(Direction face, int index, real *pos, real sigma, ITIM *itim,int reset_flag,int ** gmx_index_phase, t_topology * top, rvec * x0){
	int * p;          
	int ccc=0;
	static int flag=0;
	static int partn=0;
	real res[3];
	static int counter=0; // this counts touched mesh lines
	struct kdres *presults;
	if(reset_flag==1) { flag=0; return 0; }
	if(reset_flag==2) { counter=0; partn=0; return 0; }
	if(counter == itim->mesh.nelem)  return 0;	
	itim->alpha_index  = (int *) realloc (itim->alpha_index,(itim->n[0])* sizeof (int)); 
	itim->alphapoints  = (real *) realloc (itim->alphapoints,(itim->n[0])* sizeof (real)*3);
	itim->gmx_alpha_id = (int *) realloc (itim->gmx_alpha_id,(itim->n[0])* sizeof (int));
	partn++;
        /* find the lines within range from our particle position pos */
	presults = kd_nearest_range( itim->mesh.tree, pos, sigma + itim->alpha); 
        while( !kd_res_end( presults ) ) {
	  int test=0;
          /* get the data and position of the current result item */
          p =  (int*)kd_res_item(presults, res); /* they were init'd all to DIR_POSITIVE. Once we find a 
	  					    touching gridline, we flip the value (see below). */
	  if (*p==face) {
		if (is_in_correct_phase(index)) { // in  main cluster, layer not yet assigned OR > LAYER_OFFSET
			if(!itim->bMol)
				itim->mask[index]=itim->current_layer;
			if(itim->bMol && itim->mask[index] < LAYER_OFFSET) {  // let's add all other atoms in the molecule to the layer.
				int mol = itim->indexm[index];
			        int i1 = itim->backindex[mol];	
			        int i2 = itim->backindex[mol+1];	
				if(i2==-1) i2=itim->n[SUPPORT_PHASE]; // i2==-1 in case of last atom.
				int k ; 
	//			printf("marking atom %d, mol %d -> %d\n",index, mol,itim->current_layer);
				for(k=i1;k<i2;k++){
					if(itim->mask[k] != itim->current_layer ){
						itim->mask[k] = itim->current_layer + LAYER_OFFSET;// this way they are not tagged and will undergo surface identification.
		//				printf("marking atom %d, mol %d -> %d + OFF\n",index, mol,itim->current_layer);
					}
				}
				itim->mask[index]=itim->current_layer; // but the initial atom should be indeed tagged...
		        }
			if(itim->current_layer==1 && itim->mask[index]<LAYER_OFFSET) { 
		          itim->alpha_index[itim->nalphapoints] = index;
			  if(itim->alphapoints==NULL) exit(printf("Error reallocating alphapoints\n"));
			  if(itim->com_opt[SUPPORT_PHASE]==0) { 
	                         itim->gmx_alpha_id[itim->nalphapoints] =  gmx_index_phase[SUPPORT_PHASE][index];
		                 itim->alphapoints[3*itim->nalphapoints+0] = itim->phase[SUPPORT_PHASE][3*index+0]; 
		                 itim->alphapoints[3*itim->nalphapoints+1] = itim->phase[SUPPORT_PHASE][3*index+1];
		                 itim->alphapoints[3*itim->nalphapoints+2] = itim->phase[SUPPORT_PHASE][3*index+2];
			         itim->nalphapoints++;
			  } else {
			  	int ind,dir,k;
			  	real com[3]={0.0, 0.0, 0.0},mass=0,tot_mass=0,tmpcom[3],firstatom[3];
			  	/* in case the atoms belong to a molecule we tagged already, skip this point */
                                  for(ind=0;ind<itim->nalphapoints;ind++){ 
			  		if(itim->gmx_alpha_id[ind] == index ){
			  			ind = -1; break;
			  		}
			  	}
			  	for(k=  0 ;  k < itim->com_opt[SUPPORT_PHASE] ; k++) { 
                                      /* let's go through all atoms in the tagged molecule ...  */
			  	    ind = itim->com_opt[SUPPORT_PHASE] * ( index / itim->com_opt[SUPPORT_PHASE] ) + k;
			              /* adding them to our list ...  */		
			  	    itim->gmx_alpha_id[itim->com_opt[SUPPORT_PHASE]*itim->nalphapoints+k] =  gmx_index_phase[SUPPORT_PHASE][ ind ]; 
			  	    mass = itim->masses[gmx_index_phase[SUPPORT_PHASE][ind]];
			  	    tot_mass += mass;
                                      /* ...and computing the com. */
			  	    for(dir=0; dir<3; dir++){
			  		tmpcom[dir] = x0[itim->gmx_alpha_id[ itim->com_opt[SUPPORT_PHASE] * 
                                                        itim->nalphapoints + k ]][dir] ;
                                          if(k==0){
			  		   firstatom[dir] = tmpcom[dir];
			  	        } else { 
			  		   while(tmpcom[dir] -firstatom[dir] >  itim->box[dir]/2.){ tmpcom[dir] -= itim->box[dir];} 
			  		   while(tmpcom[dir] -firstatom[dir] < -itim->box[dir]/2.){ tmpcom[dir] += itim->box[dir];}
                                          }
			  		com[dir] += tmpcom[dir] * mass;
                                      }
			  	}
			  	
			  	for(dir=0;dir <3; dir++){
			  	        com[dir]/=tot_mass;
			  		while(com[dir] >  itim->box[dir]/2.) com[dir] -= itim->box[dir];
			  		while(com[dir] < -itim->box[dir]/2.) com[dir] += itim->box[dir];
			  	        itim->alphapoints[3*itim->nalphapoints+dir]=com[dir];
			  	}
			          itim->nalphapoints++;
			  }
			}
		}	
		*p=-1*face;  // flip  the flag of the testline.
		counter ++;
		if(counter == itim->mesh.nelem) {
                       /* Now all the testline have been associated to molecules, and have therefore all been flipped */
		      kd_res_free( presults );
                      counter = 0 ;  return 0 ;
                }
	  }
          /* go to the next entry */
          kd_res_next( presults );
        }
	kd_res_free( presults );
	return counter;
}

int check_itim_testlines_periodic(Direction face, int index, real * pos,real sigma, ITIM *itim, int ** gmx_index_phase,t_topology * top, rvec * x0){
	int ret=1;
	real periodic[2];
	real border_positive[2];
	real border_negative[2];
	border_positive[0]=itim->box[0]/2. - sigma - itim->alpha ;
	border_positive[1]=itim->box[1]/2. - sigma - itim->alpha ;
	border_negative[0]=-itim->box[0]/2. + sigma + itim->alpha ;
	border_negative[1]=-itim->box[1]/2. + sigma + itim->alpha ;
	ret = check_itim_testlines(face,index, pos,sigma,itim,0,gmx_index_phase,top, x0);
	if(ret==0) return 0;

	if(pos[0] <= border_negative[0] || pos[0] >= border_positive[0] ){
		if(pos[0]<= border_negative[0]){
			periodic[0]=pos[0]+itim->box[0];
		} else { 
			periodic[0]=pos[0]-itim->box[0];
		}
		periodic[1]=pos[1];
		ret = check_itim_testlines(face,index,periodic,sigma,itim,0,gmx_index_phase,top,x0);
		if(ret==0) return 0;
	}
	if(pos[1] <= border_negative[1] || pos[1] >= border_positive[1] ) { 
		if(pos[1]<= border_negative[1]){
			periodic[1]=pos[1]+itim->box[1];
		} else { 
			periodic[1]=pos[1]-itim->box[1];
		}
		periodic[0]=pos[0];
		ret = check_itim_testlines(face,index, periodic,sigma,itim,0,gmx_index_phase,top,x0);
		if(pos[0] <= border_negative[0] || pos[0] >= border_positive[0] ){
			if(pos[0]<= border_negative[0]){
				periodic[0]=pos[0]+itim->box[0];
			}
			else { 
				periodic[0]=pos[0]-itim->box[0];
			}
			ret = check_itim_testlines(face,index, periodic,sigma,itim,0,gmx_index_phase,top,x0);
			if(ret==0) return 0;
		}
        }
	check_itim_testlines((Direction)0,0, NULL,0,NULL,1,NULL,NULL,NULL); // this just resets the flag within check_itim_testlines.
	return ret;
}

int projection_positive(const void * a, const void * b) {
	real a2,b2; //NOTE this is weighted
	int i,j;
	ITIM * itim = global_itim;
	i = *(int*)a;	
	j = *(int*)b;	
	a2 = itim->phase[SUPPORT_PHASE][3*i+2] + itim->radii[i] ;
	b2 = itim->phase[SUPPORT_PHASE][3*j+2] + itim->radii[j] ;
	   // here we use groups 1 & 2 as opposite phases to determine the inclusive set.
	   // mask just stacks all groups masks, one after each other, so we can use this simple
	   // trick here. The support group in this case has to coincide with groups 1+2, with no overlaps.
	   // TODO: Add a test about it.
        if(! is_in_correct_phase(i)) a2-=10000.; // not in the biggest cluster, let's shift it back away
        if(! is_in_correct_phase(j)) b2-=10000.; 
	return (a2>b2?-1:1);
}
int projection_negative(const void * a, const void * b) {
	real a2,b2; //NOTE this is weighted
	int i,j;
	ITIM * itim = global_itim;
	i = *(int*)a;	
	j = *(int*)b;	
	a2 = itim->phase[SUPPORT_PHASE][3*i+2] - itim->radii[i] ;
	b2 = itim->phase[SUPPORT_PHASE][3*j+2] - itim->radii[j] ;
        if(! is_in_correct_phase(i)) a2+=10000.; // not in the biggest cluster, let's shift it back away
        if(! is_in_correct_phase(j)) b2+=10000.; 
	return (a2<b2?-1:1);
}

int projection(const void * a, const void * b) {
	real a2,b2; //NOTE this is weighted
	int i,j,side;
	ITIM * itim = global_itim;
	side=itim->side;
	i = *(int*)a;	
	j = *(int*)b;	
	if(side*side!=1) exit(printf("Internal error, side must me either 1 or -1. %s:%d",__FILE__,__LINE__));
	a2 = itim->phase[SUPPORT_PHASE][3*i+2] + itim->radii[i] ;
	b2 = itim->phase[SUPPORT_PHASE][3*j+2] + itim->radii[j] ;
	if(itim->bInclusive){
	   int size[3];
	   size[0]=itim->n[SUPPORT_PHASE];
	   size[1]=itim->n[SUPPORT_PHASE+1];
	   size[2]=itim->n[SUPPORT_PHASE+2];
	   // here we use groups 1 & 2 as opposite phases to determine the inclusive set.
	   // mask just stacks all groups masks, one after each other, so we can use this simple
	   // trick here. The support group in this case has to coincide with groups 1+2, with no overlaps.
	   // TODO: Add a test about it.
	   if(i<size[1]){// "main group", of which we are taking atoms from the biggest cluster
               if(itim->mask[i+size[0]]==-1||itim->mask[i+size[0]]>0) a2-=10000.; // not in the biggest cluster, let's shift it back away
           } else { 
               if(itim->mask[i+size[0]]==0||itim->mask[i+size[0]]>0) a2-=10000.; // in the biggest cluster
	   }
	   if(j<size[1]){// same for atom j
               if(itim->mask[j+size[0]]==-1||itim->mask[j+size[0]]>0) b2-=10000.;
           } else { 
               if(itim->mask[j+size[0]]==0||itim->mask[j+size[0]]>0) b2-=10000.;
           }
	} else { 
           if(itim->mask[i]==-1||itim->mask[i]>0) a2-=side*10000.;
           if(itim->mask[j]==-1||itim->mask[j]>0) b2-=side*10000.;
        }
	return (a2>b2?-side:side);
}

// TODO !! NOTE: this allows only computation of surface molecules of the SUPPORT_PHASE
int compute_itim_points(Direction direction, ITIM *itim, int ** gmx_index_phase,t_topology * top, rvec * x0,int * mask){
	int i=0,index,result;
#ifdef TIME_PROFILE
        struct timeval tp;
        struct timeval tp2;
        gettimeofday(&tp, NULL);
#endif
	// sorting is done only on the index, actual positions are left untouched.
	itim->side=1;
	qsort((void*)itim->phase_index[SUPPORT_PHASE], itim->n[SUPPORT_PHASE], sizeof(int) , projection_positive);
#ifdef TIME_PROFILE
        gettimeofday(&tp2, NULL);
        fprintf(stdout,"Time to quicksort for itim: millisec=%f\n",1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
        /* The positive side */
        i=0;
        result=1;
	do { 
		index = itim->phase_index[SUPPORT_PHASE][i];
		if(is_in_correct_phase(index)){
			result = check_itim_testlines_periodic(DIR_POSITIVE,index,&itim->phase[SUPPORT_PHASE][3*index],itim->radii[index],itim,gmx_index_phase,top,x0) ; 
		}
          	i++ ; 
		if(i>=itim->n[SUPPORT_PHASE]) {printf("Error: all (%d) particles scanned, but did not associate all testlines on the positive side...\n",itim->n[SUPPORT_PHASE]); return 0; }
	} while (result) ;
	check_itim_testlines((Direction)0,0, NULL,0,NULL,1,NULL,NULL,NULL); 
	check_itim_testlines((Direction)0,0, NULL,0,NULL,2,NULL,NULL,NULL); 
        i = 0;
        result=1;
	itim->side=-1;
	qsort((void*)itim->phase_index[SUPPORT_PHASE], itim->n[SUPPORT_PHASE], sizeof(int) , projection_negative);
	do { 
		index = itim->phase_index[SUPPORT_PHASE][i];
		if(is_in_correct_phase(index)){
			result = check_itim_testlines_periodic(DIR_NEGATIVE,index,&itim->phase[SUPPORT_PHASE][3*index],itim->radii[index],itim,gmx_index_phase,top,x0) ; 
		}
          	i++; 
		if(i>=itim->n[SUPPORT_PHASE]) {printf("Error: all (%d) particles scanned, but did not associate all testlines on the negative side...\n",itim->n[SUPPORT_PHASE]); return 0; }
	} while (result) ;
	check_itim_testlines((Direction)0,0, NULL,0,NULL,1,NULL,NULL,NULL); 
	check_itim_testlines((Direction)0,0, NULL,0,NULL,2,NULL,NULL,NULL); 
	for (int i=0 ; i < itim->n[SUPPORT_PHASE];i++) itim->mask[i] %= LAYER_OFFSET;
	return 1;
}


void arrange_alpha_points(ITIM *itim, int ** gmx_index_phase, t_topology * top,rvec * x0){
	int i ;

	if(Surface!=NULL) {  kd_free(Surface); Surface=NULL; }
	switch(itim->geometry) {
		case SURFACE_PLANE:  Surface  = kd_create( 2 ); break;
		default: exit(printf("Geometry not implemented [arrange_alpha_points()]\n"));
	}
	switch(itim->method){
	      case METHOD_ITIM : 
                    /* Note: for ITIM the coms have been already computed in the check_itim_testlines() function */
		    for(i=0;i<itim->nalphapoints;i++) 
			  kd_insert(Surface, &itim->alphapoints[3*i], &itim->alphapoints[3*i+2]);
	      break;
              default: exit(printf("Method not implemented\n"));
         }
}

real D2dist (real *a, real *b){
	return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
}



int distance(void* ref, const void * a, const void * b) {
	real d1,d2;
	real a0,a1,b0,b1,ref0,ref1;
	a0=((real*)a)[0]; a1=((real*)a)[1]; b0=((real*)b)[0]; b1=((real*)b)[1];
	ref0=((real*)ref)[0]; ref1=((real*)ref)[1]; 
	d1= (a0-ref0)*(a0-ref0) +
	    (a1-ref1)*(a1-ref1) ;
	d2= (b0-ref0)*(b0-ref0) +
	    (b1-ref1)*(b1-ref1) ;
        return  ((d1>d2)?1:-1);
}

int periodic_distance(const void * a, const void * b) {
	real d1,d2;
	real a0,a1,b0,b1,ref0,ref1;
	real box[3];
	box[0]=global_itim->box[0];
	box[1]=global_itim->box[1];
	box[2]=global_itim->box[2];
	ref0  =global_real_pointer[0];
	ref1  =global_real_pointer[1];
	a0=((real*)a)[0]; a1=((real*)a)[1]; b0=((real*)b)[0]; b1=((real*)b)[1];

	a0 -= ref0; a1 -= ref1; b0 -= ref0; b1 -= ref1;
	while (a0> box[0]/2.) a0 -= box[0]; while (a0<-box[0]/2.) a0 += box[0];
	while (a1> box[1]/2.) a1 -= box[1]; while (a1<-box[1]/2.) a1 += box[1];
	while (b0> box[0]/2.) b0 -= box[0]; while (b0<-box[0]/2.) b0 += box[0];
	while (b1> box[1]/2.) b1 -= box[1]; while (b1<-box[1]/2.) b1 += box[1];
	d1= a0*a0 + a1*a1 ;
	d2= b0*b0 + b1*b1 ;
        return  ((d1>d2)?1:-1);
}

void dump_phase_points(PHASE phase, t_topology *top, FILE* cid){
	int i,atom_index=0;
	char phase_name[16];
	ITIM * itim=global_itim;
	switch(phase){
		case INNER_PHASE: sprintf(phase_name,"inner"); break;
		case OUTER_PHASE: sprintf(phase_name,"outer"); break;
		default: exit(printf("Internal error in dump_phase_points()\n"));
	}
        fprintf(cid,"%s\n",*(top->name));
        fprintf(cid,"%d\n",itim->n[phase]);
 	for(i=0;i<itim->n[phase];i++){
                atom_index=itim->gmx_index[phase][i];
                fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                            top->atoms.atom[atom_index].resind+1,
                            *(top->atoms.resinfo[top->atoms.atom[atom_index].resind].name),
                            *(top->atoms.atomname[atom_index]),
			    atom_index+1,
			    itim->phase[phase][3*i],
			    itim->phase[phase][3*i+1],
			    itim->phase[phase][3*i+2],
                            0.0,0.0,0.0);
	}
        fprintf(cid,"%f %f %f\n",itim->box[0],itim->box[1],itim->box[2]);
}

void spol_atom2molindex(int *n, int *index, int*backindex, t_block *mols)
{
    int nmol, i, j, m;

    nmol = 0;
    i    = 0;
    m=0;
    while (m < mols->nr) {
	backindex[m]=-1;
	m++;
    }

    m = 0;
    for (i =0 ; i < *n ; i++) { 
        while (m < mols->nr &&  ! ( mols->index[m] <= index[i]  && index[i] < mols->index[m+1] ) )  { 
		m++;
	}
        if (m == mols->nr)
        {
            gmx_fatal(FARGS, "index[%d]=%d does not correspond to the first atom of a molecule", i+1, index[i]+1);
        }
	if(backindex[m]==-1) backindex[m]=i; 
        index[i]=m;	
    }
}

void compute_normal(rvec p1,rvec p2,rvec p3,real * normal,GEOMETRY geometry){
            rvec v1,v2,vn;
            rvec_sub(p2,p1,v1);                                  
            rvec_sub(p3,p1,v2);                                  
            cprod(v1,v2,vn);
            unitv(vn,vn);
            /* SAW: here we are assuming that we are using the -center option,
                    which should anyway always on for the planar case and spherical cases... Fix this. */
            if (geometry==SURFACE_PLANE)
                 if( vn[2]*p1[2]<0 ) vn[2]*=-1;
            normal[0]=vn[0]; normal[1]=vn[1]; normal[2]=vn[2];
}

real perform_interpolation( struct kdtree *Surface, real * P, real * normal, ITIM * itim, int i, int phase) {
		real p1[3],p2[3],p3[3],diff[3];
                rvec v1,v2,vn;
		real a1,a2,a3,atot,dist;
                struct kdres * presults;
		int j=0,found=0;
		real * zpos;	

		/* p1 */
		presults = kd_nearest_range( Surface, P,0.8+global_itim->alpha); 
		if(kd_res_end(presults)){  printf("Did not find any neighbor on the surface, skipping...\n"); return 99.999;}
		zpos    = (real*) kd_res_item(presults, p1); 
		if(kd_res_end(presults)) exit(printf("Error, no projected surface point found\n"));
		p1[2]=*zpos;
		while(p1[2]*P[2]<0){ /*i.e. they are on opposite sides, this can happen only with ITIM*/
			kd_res_next( presults );
			if(kd_res_end(presults)){ if(NULL!=getenv("GITIM_WARNING")) fprintf(stdout,"Warning: no suitable point for interpolation found \n"); return 99.999 ; }
			zpos    = (real*) kd_res_item(presults, p1); 
			p1[2]=*zpos;
		}

                rvec_sub(p1,P,diff); if(norm2(diff)<1e-6) { kd_res_free(presults); return 0.0 ;  } 
                do {
			kd_res_next( presults );
			if(kd_res_end(presults)){ if(NULL!=getenv("GITIM_WARNING")) fprintf(stdout,"Warning: no suitable point for interpolation found \n"); return 99.999 ; }
                	kd_res_item(presults, p2); 

			zpos    = (real*) kd_res_item(presults, p2); 
			p2[2]=*zpos;

		} while(p2[2]*P[2]<0);

                rvec_sub(p2,P,diff); if(norm2(diff)<1e-6) { kd_res_free(presults); return 0.0 ;  } 

		/*Now start iterating over all other (sorted) surface atoms p3, to find a 
		  triangle which p1-p2-p3 encloses our point P, i.e. s.t. a(123)=a(124)+a(234)+a(134)*/
		a1=fabs(triangle_area(p1,p2,P,itim->box));
		while (kd_res_next( presults )) {  
		    zpos    = (real*) kd_res_item(presults, p3); 
	            p3[2]=*zpos;

		    if(p3[2]*P[2]<0) continue; /* they are not on the same side */
		    /* triangle_area returns the area of the triangle, projected on the xy plane */
		    atot=fabs(triangle_area(p1,p2,p3,itim->box));
		    a2=fabs(triangle_area(p1,p3,P,itim->box));
		    a3=fabs(triangle_area(p2,p3,P,itim->box));
		    if(fabs(a1+a2+a3 - atot) < 1e-6 ){

				dist = interpolate_distance(p1,p2,p3,P); /* this is not really a distance, as it is signed, 
									     but we carry the sign to have the complete information, 
                                                                             as particles located behind the surface might exist */
				dist = P[2]-dist; 
				if(P[2]<0.) dist = -dist ; /* the other side...*/

		     		compute_normal(p1,p2,p3,normal,itim->geometry);
                                kd_res_free(presults);
				return  dist;
			} 
                }
			
                /* no suitable triangle found */
                normal[0]=normal[1]=normal[2]=NORMAL_UNDEFINED;
	        dist = P[2] -p1[2];
		if(P[2]<0.) dist = -dist ; /* the other side...*/

                kd_res_free(presults);
		return dist;
}

void dump_slabs(t_topology* top, int dump_mol){
	int i,atom_index=0, phase_index, layer;
	static int frame[2]={0,0};
	static FILE * (fp[2])={NULL,NULL};
        ITIM * itim=global_itim; 
	int isizem, *indexm=NULL,*backindex=NULL;
	char fname[2048];
	sprintf(fname,"layers.pdb");
	if(dump_mol)sprintf(fname,"layers_mol.pdb");
	frame[dump_mol]++;

	if(fp[dump_mol]==NULL){ 
		fp[dump_mol] = fopen(fname,"w"); 
	} else { 
		fclose(fp[dump_mol]);
		fp[dump_mol] = fopen(fname,"a"); 
	} 
	fprintf(fp[dump_mol],"REMARK    GENERATED BY A COMPUTER\n");
	fprintf(fp[dump_mol],"TITLE     LAYERS t= %d\n",frame[dump_mol]);
	fprintf(fp[dump_mol],"REMARK    THIS IS A SIMULATION BOX\n");
	fprintf(fp[dump_mol],"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%10s%3d\n",
			  itim->box[0]*10.,itim->box[1]*10.,itim->box[2]*10.,
			   90.,90.,90.,"P  1 ",1);
	fprintf(fp[dump_mol],"MODEL         %d\n",frame[dump_mol]);

	if(dump_mol) { 
		isizem=itim->n[SUPPORT_PHASE];
		snew(indexm,isizem);
		snew(backindex,top->mols.nr+1);
		for(i=0;i<itim->n[SUPPORT_PHASE];i++) indexm[i]=itim->gmx_index[SUPPORT_PHASE][i]; // NOTE: TODO check "additional"  under the PATCH case
                spol_atom2molindex(&isizem, indexm,backindex, &(top->mols));
        }
 	for(i=0;i<itim->n[SUPPORT_PHASE];i++){
	       		
	        int mol,i1,i2;
		real dist;
	        phase_index= i ; 
		real normal[3];
	        atom_index = itim->gmx_index[SUPPORT_PHASE][phase_index];

		dist = perform_interpolation(Surface, &(itim->phase[SUPPORT_PHASE][3*i]), &(normal[0]),itim,i,SUPPORT_PHASE);
//	        if (dump_mol){
//			int minlayer = 100000;
//			int k;
//	        	mol  = indexm[i];
//	        	i1 = backindex[mol] ;
//	        	i2 = backindex[mol+1];
//			if(i2==-1) i2=itim->n[SUPPORT_PHASE]; // i2==-1 in case of last atom.
//			for(k=i1;k<i2;k++){
//				if (itim->mask[k] >0 && itim->mask[k]%LAYER_OFFSET < minlayer ) { 
//					minlayer = itim->mask[k] % LAYER_OFFSET;
//				}
//			}	
//			layer= ( minlayer==100000 ? 0 :minlayer) ;	
		//	printf("mol %d -> %d %d %d %d\n",mol,i1,i2,layer,itim->mask[i]);
	//	} else {
			layer = itim->mask[i] % LAYER_OFFSET;
 	//	}
	        if(layer< itim->maxlayers+1){
                  fprintf(fp[dump_mol],"%-6s%5d %4s%1s%3s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
	         	   "ATOM",
	         	   atom_index+1,
                            *(top->atoms.atomname[atom_index]),
	         	   " ",
                            *(top->atoms.resinfo[top->atoms.atom[atom_index].resind].name),
	         	   " ",
	         	   (top->atoms.atom[atom_index].resind+1)%10000, // pdb accepts at most 4-digit residue numbers
	         	   " ",
	        	   10*itim->phase[SUPPORT_PHASE][3*i],
	        	   10*itim->phase[SUPPORT_PHASE][3*i+1],
	        	   10*itim->phase[SUPPORT_PHASE][3*i+2],
			   (double)layer,
                            10.*dist,
	         	   " ",
	         	   itim->phase[SUPPORT_PHASE][3*i+2]>0?"1":"2");
                }
        }
        fprintf(fp[dump_mol],"TER\nENDMDL\n");

	if(dump_mol) { 
		free(indexm);
		free(backindex);
        }
}

void dump_surface_points(t_topology* top,FILE* cid){
	int i,atom_index=0;
        ITIM * itim=global_itim; 
        fprintf(cid,"%s\n",*(top->name));
        fprintf(cid,"%d\n",itim->nalphapoints);
 	for(i=0;i<itim->nalphapoints;i++){
               atom_index = itim->gmx_alpha_id[i];
// printf("%s ----> %d %d\n",*(top->atoms.atomname[atom_index]),atom_index,i);
               fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                           top->atoms.atom[atom_index].resind+1,
                           *(top->atoms.resinfo[top->atoms.atom[atom_index].resind].name),
                           *(top->atoms.atomname[atom_index]),
	       	    atom_index+1,
	       	    itim->alphapoints[3*i],
	       	    itim->alphapoints[3*i+1],
	       	    itim->alphapoints[3*i+2],
                           0.0,0.0,0.0);
        }
        fprintf(cid,"%f %f %f\n",itim->box[0],itim->box[1],itim->box[2]);
}
/* TODO: this works only if the SUPPORT and INNER phases are the same (in the input). FIX it by allowing to compute the molecules corresponting to the SUPPORT phase. This needs propoer indexing */
void dump_surface_molecules(t_topology* top,FILE* cid,atom_id ** gmx_index_phase){
	int i,j,atom_index=0,residue_index=0,surface_index=0;
	static int resNR=0;
	static int atomsinRes=0;

        int *residuelist;
	int checkedresidues=0;
        ITIM * itim=global_itim; 

        residuelist=(int*)malloc(itim->n[SUPPORT_PHASE]);

        atomsinRes=itim->dump_mol;
 	for(i=0;i<itim->nalphapoints;i++){
               surface_index = itim->alpha_index[i];
               int offset = surface_index%atomsinRes;//(itim->gmx_alpha_id[i]%atomsinRes);
               surface_index -= offset;  /*  here we rewind to the beggining of the residue ... */
               atom_index = itim->gmx_index[SUPPORT_PHASE][surface_index];
               for(j=0;j<checkedresidues;j++){
		   if(residuelist[j]==atom_index) break;
               }
               if(j==checkedresidues) { // otherwise, we have already printed this molecule
                 residuelist[checkedresidues]=atom_index; 
                 checkedresidues++;
               }
        }

        fprintf(cid,"%s\n",*(top->name));
       fprintf(cid,"%d\n",checkedresidues*atomsinRes); 
        checkedresidues=0;
							   
 	for(i=0;i<itim->nalphapoints;i++){
               surface_index = itim->alpha_index[i];
               int offset = surface_index%atomsinRes;//(itim->gmx_alpha_id[i]%atomsinRes);
 //              if(offset!=0) printf("Here:%d %d %d\n",i,offset,surface_index);
               surface_index -= offset;  /*  here we rewind to the beggining of the residue ... */
               atom_index = gmx_index_phase[SUPPORT_PHASE][surface_index];
               for(j=0;j<checkedresidues;j++){
		   if(residuelist[j]==atom_index) break;
               }
               if(j==checkedresidues) { // otherwise, we have already printed this molecule
                 residuelist[checkedresidues]=atom_index; 
                 checkedresidues++;
 //printf("%s ----> %d %d %d\n",*(top->atoms.atomname[atom_index]),atom_index,checkedresidues,i);
	         for(j=0; j< atomsinRes ; j++) {
                    fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", /* if you gdb up till here, remember that 
									           this broken part of code requires 
									           INNER_PHASE to be the same as SUPPORT_PHASE */
                           top->atoms.atom[atom_index+j].resind+1,
                           *(top->atoms.resinfo[top->atoms.atom[atom_index + j].resind].name),
                            *(top->atoms.atomname[atom_index+j]),
	       	    atom_index+j+1,
	       	    itim->phase[SUPPORT_PHASE][3*(surface_index+j)+0],
	       	    itim->phase[SUPPORT_PHASE][3*(surface_index+j)+1],
	       	    itim->phase[SUPPORT_PHASE][3*(surface_index+j)+2],
                           0.0,0.0,0.0);
                 }             
              }
        }
        fprintf(cid,"%f %f %f\n",itim->box[0],itim->box[1],itim->box[2]);
        free(residuelist);
}




void init_itim(int nphases,int nadd_index) { 
	ITIM * itim;
	int i;
	global_itim=(ITIM*)malloc(sizeof(ITIM));
	itim=global_itim;
	sprintf(itim->method_name[METHOD_ITIM],"itim");
	sprintf(itim->method_name[METHOD_A_SHAPE],"alpha-shapes");
        itim->masses=global_masses;
        itim->charges=global_charges; 
	if(itim->masses==NULL) exit(printf("Internal error, itim->masses not allocated\n"));
	itim->mesh.nelem=0;
	itim->mesh.tree=NULL;
	itim->mesh.flag=NULL;
	itim->nphases = nphases;
	itim->RANDOM_PHASE= nphases-1;
        itim->current_layer = 1;
        itim->ngmxphases=itim->nphases-1;
	itim->side=1;

	itim->n = (int*)malloc(nphases*sizeof(int));
	itim->mask_add = (int**)malloc(nadd_index*sizeof(int*));
        itim->phase_index =  (int**)malloc(nphases*sizeof(int*));
        itim->inclusive_map =  (int*)malloc(2*nphases*sizeof(int));
        itim->phase =  (real **)malloc(nphases*sizeof(real*));
        itim->periodic =  (PERIODIC*)malloc(nphases*sizeof(PERIODIC));

	itim->nalphapoints = 0;
	itim->alphapoints = NULL;
        itim->alpha_index= NULL; 
        itim->radii = NULL;
        itim->gmx_alpha_id= NULL;
	for(i=0;i<nphases;i++){
        	itim->n[i] = 0;
	        itim->phase_index[i]= NULL ;
	        itim->phase[i] = NULL ;
		itim->periodic[i]=NONE; 
	}
	for(i=0;i<nadd_index;i++){
		itim->mask_add[i]=NULL;
	}
	for(i=0;i<64;i++){
                itim->com_opt[i]=0;
        }
        itim->dump_mol=0;
	itim->periodic[INNER_PHASE]=NONE;
	itim->dump_surface_points=dump_surface_points;
	itim->dump_surface_molecules=dump_surface_molecules;
	itim->dump_phase_points=dump_phase_points;
        
}

void generate_mask_ns(int bCluster, int bInclusive, rvec * gmx_coords, int ** mask, int *nindex, int ** gmx_index_phase, matrix box, int  ng, int natoms, t_pbc pbc){
/* New mask: -1 -> not in the main cluster
	      0 -> layer not assigned
              i -> i-th layer (positive or negative) */

    	int        ** cluster_size;
    	int        ** cluster_map;
    	int        ** cluster_analyzed;
    	int        ** cluster_index;
	int nclusters=-1;
	int   totsize=0, * incsize;

	incsize = (int*) malloc(sizeof(int)*ng);
	for(int g=0;g<ng;g++){
		incsize[g]=totsize;
		totsize+=nindex[g];
	}
	if(*mask==NULL){
	  *mask=(int*)malloc(totsize*sizeof(int));
        }
	if(!bCluster){
	      // by default, use all atoms;
             for(int atom=0 ; atom <  totsize ; atom++) (*mask)[atom]=0;
	     free(incsize);
             return ;
        } else {
             for(int atom=0 ; atom <  totsize ; atom++) (*mask)[atom]=-1;
        }


    	snew(cluster_size,  	ng);
    	snew(cluster_map,   	ng);
    	snew(cluster_index, 	ng);
    	snew(cluster_analyzed,  ng);
        for(int phase=SUPPORT_PHASE; phase < ng; phase++) { 
    		snew(cluster_size[phase], nindex[phase]);
    		snew(cluster_map[phase], nindex[phase]);
    		snew(cluster_analyzed[phase],nindex[phase]);
    		snew(cluster_index[phase],nindex[phase]);
	}
	totsize=0;
        for (int g = 0; g < ng; g++)
        {
		if(bInclusive && g==0 ) continue; /* we don't need it, we'll copy the info in reinit_mask() */ 
		struct kdtree *kd =  kd_create(3);
		struct kdres *kdresults;
		rvec  xx;
		real *idata;
		idata = (real*) malloc(sizeof(real)*nindex[g]);
		for(int i=0;i<nindex[g];i++){
		    idata[i]=i;
		    copy_rvec(gmx_coords[gmx_index_phase[g][i]],xx);
		    kd_insert3p(kd,xx[0],xx[1],xx[2],&idata[i],CLUSTERCUT[g]);
		    //kd_insert3(kd,xx[0],xx[1],xx[2],&idata[i]);
		}
	        for (int i = 0; i < nindex[g]; i++){
			cluster_size[g][i]=0;
			cluster_map[g][i]=i;
	                cluster_analyzed[g][i] = 0; 
		}
		nclusters=0;
	        int j=0,ii,ij;
	        int nanalyzed=1;
		rvec pos;
		real coord[3];
	        for (int i = 0; i < nindex[g]; i++)
	        {
	            if (cluster_analyzed[g][i] == 0) {
	             cluster_analyzed[g][i]=1;
		     cluster_map[g][i]=nclusters;
	             cluster_size[g][nclusters]++; 
	             cluster_index[g][j]=i;
	             j++;
	             while (j>=nanalyzed){ 
	                  ii = cluster_index[g][nanalyzed-1];
	                  //gmx_ana_nbsearch_first_within(data,gmx_coords[gmx_index_phase[g][ii]],&ij);
			  coord[0]=gmx_coords[gmx_index_phase[g][ii]][0];
			  coord[1]=gmx_coords[gmx_index_phase[g][ii]][1];
			  coord[2]=gmx_coords[gmx_index_phase[g][ii]][2];
			  kdresults = kd_nearest_range( kd, coord, CLUSTERCUT[g]); 
	                  while(!kd_res_end( kdresults )){
	// put a type-based cutoff analysis
				ij = (int)(*((real*)kd_res_item( kdresults, pos )));
{
				real x2[3];
				real d=0.0;
				for(int i =0 ; i < 3 ; i++) { 
					real delta;
					x2[i]=gmx_coords[gmx_index_phase[g][ij]][i];
					delta = coord[i]-x2[i];
					while(delta>global_itim->box[i]/2.) delta-=global_itim->box[i];
					while(delta<-global_itim->box[i]/2.) delta+=global_itim->box[i];
					d+=delta*delta;
				}
}
	                     	if (cluster_analyzed[g][ij] == 0 ) { 
				       //printf("g=%d j=%d nan=%d\n",g,j,nanalyzed);
				       cluster_map[g][ij]=nclusters;
				       cluster_size[g][nclusters]++;
	                               cluster_analyzed[g][ij]=1;
				       //printf("cluster_index[%d][%d]=%d\n",g,j,cluster_index[g][j]);
	                               cluster_index[g][j]=ij;
	                               j++;
				}
				kd_res_next( kdresults );
			           //gmx_ana_nbsearch_next_within(data,&ij);
	                  }
		 	  kd_res_free( kdresults );
	                  nanalyzed++;
	             }
		     nclusters++; 
	            } 
	         }

		 {
		   int msize=0;
		   int largest=-1;
		   for(int i = 0 ; i< nindex[g] ; i++) if(cluster_size[g][i]>=msize){msize=cluster_size[g][i] ; largest = i;}
		   for(int i = 0 ; i< nindex[g] ; i++) if(cluster_map[g][i]==largest){(*mask)[incsize[g]+i]=0 ; }
		   if(msize<=1){
			 if(getenv("NO_CLUSTERSIZE_ERROR") == NULL)
			   exit(printf("ERROR: the largest cluster for group index %d is %d, \nyou might have choosen a too small cutoff for the cluster algorithm\nSet the variable NO_CLUSTERSIZE_ERROR=1 to suppress this\n",g,msize));
                   }
		   free(idata);
		// gmx_ana_nbsearch_free(data);
                   kd_free( kd );
		 }
	}
        for(int phase=SUPPORT_PHASE; phase < ng; phase++) { 
    		free(cluster_size[phase]);
    		free(cluster_map[phase]);
    		free(cluster_analyzed[phase]);
    		free(cluster_index[phase]);
	}
    	free(cluster_size);
    	free(cluster_map);
    	free(cluster_index);
    	free(cluster_analyzed);
}





void reinit_mask(ITIM * itim){
	if(itim->bInclusive){
	   int atom;
	   int offset=0,phase;
           for(phase=INNER_PHASE;phase <=OUTER_PHASE ; phase++) { 
		offset+=itim->n[phase-1];
		for(atom=0;atom<itim->n[phase];atom++){
			/* copy from the group "phase" to the corresponding portion of the SUPPORT_GROUP map" */
			itim->mask[atom+itim->inclusive_map[2*phase]] = itim->mask[offset+atom];
		}	
           } 
	}
}

void arrange_datapoints( ITIM *itim, rvec * gmx_coords, int *nindex, int ** gmx_index_phase,int * mask){
/* NOTE:  Here we assume that particles are within a [-box/2:box/2]^3. */
	int atom,additional=0,sign,sign2;
	const real scale = 3.5; /* This defines the thickness of the periodic border, in units of alpha*/
	real radius=0;
	int realloc_factor=2;
	real x,y,z;
	int phase;
        for(phase=(int)SUPPORT_PHASE; phase < itim->ngmxphases ; phase++) { 
          /* at the first iteration, itim->n are all zero  (initialized by init_itim() )*/

          if(itim->n[phase]<nindex[phase]) { 
                itim->n[phase]=nindex[phase];
                itim->phase[phase] = (real * ) realloc(itim->phase[phase],itim->n[phase]*sizeof(real) * 3);
                itim->phase_index[phase] = (int * ) realloc(itim->phase_index[phase],itim->n[phase]*sizeof(int) );
          }
	  if(phase==SUPPORT_PHASE){
                for(int phase_add=0; phase_add < itim->nadd_index; phase_add++) { 
			if(itim->mask_add[phase_add]==NULL){
			 	itim->mask_add[phase_add] = (int*) realloc(itim->mask_add[phase_add],itim->n[SUPPORT_PHASE]*sizeof(int));
				for(int j=0;j<itim->n[SUPPORT_PHASE];j++){
					itim->mask_add[phase_add][j]=0;
					for(int k=0;k<nindex[itim->ngmxphases+phase_add];k++){
						//printf(" gmx_index_phase[%d][%d] = %d == %d = itim->phase_index[SP][%d] \n",itim->ngmxphases+phase_add,k,
//gmx_index_phase[itim->ngmxphases+phase_add][k] , gmx_index_phase[SUPPORT_PHASE][j],j);	
						if (gmx_index_phase[itim->ngmxphases+phase_add][k] == gmx_index_phase[SUPPORT_PHASE][j]){
							itim->mask_add[phase_add][j]=1;		
							break;
						} 
					}
					//printf("map[%d] (%d) = %d\n",j,gmx_index_phase[SUPPORT_PHASE][j],itim->mask_add[phase_add][j]);
				}
			}
		}
	  }
          switch (itim->periodic[phase]) { 
		case NONE:
		case FULL:
          	  switch (itim->geometry) { 
                    case SURFACE_PLANE:
                 	for(atom=0 ; atom < nindex[phase] ; atom++) { 
                 	   itim->phase[phase][3*atom] = 
			      (gmx_coords)[gmx_index_phase[phase][atom]][(itim->normal+1)%3];  
                 	   itim->phase[phase][3*atom+1] = 
			      (gmx_coords)[gmx_index_phase[phase][atom]][(itim->normal+2)%3];  
                 	   itim->phase[phase][3*atom+2] = 
				   (gmx_coords)[gmx_index_phase[phase][atom]][itim->normal];  
			
                 	   itim->phase_index[phase][atom] =  atom; 
                                /* now points[2] is the coordinate along the normal, and the other are such 
                                   that they form a positively ordered triplet (eg. xyz, yzx or zxy) */
			   
               		}
			if(itim->bInclusive){
				if(( phase==INNER_PHASE||phase==OUTER_PHASE)) {
					int support_atom;
					for(support_atom=0;support_atom<nindex[SUPPORT_PHASE]; support_atom++){
						if (gmx_index_phase[phase][0] == gmx_index_phase[SUPPORT_PHASE][support_atom] ) {
							itim->inclusive_map[2*phase]= support_atom;
						}
						if ( gmx_index_phase[phase][nindex[phase]-1] == gmx_index_phase[SUPPORT_PHASE][support_atom] ) {
							itim->inclusive_map[2*phase+1]= support_atom;
						}
					}
				}
			   }
		        break;
                    default: exit(printf("Geometry not implemented [arrange_datapoints()]\n"));
                  }
		break; // FULL || NONE
	 	case PATCH:
          	  switch (itim->geometry) { 
			case SURFACE_PLANE:
		        /*Let's provisionally allocate more space*/	
                        itim->phase[phase] = (real * ) realloc(itim->phase[phase],realloc_factor*nindex[phase]*sizeof(real) * 3);
                        itim->phase_index[phase] = (int * ) realloc(itim->phase_index[phase],realloc_factor*nindex[phase]*sizeof(int) );
            		if(phase==SUPPORT_PHASE) itim->radii=  (real *)realloc(itim->radii, realloc_factor*nindex[phase]*sizeof(real));

                 	for(atom=0 ; atom < nindex[phase] ; atom++) { 
			   x = (gmx_coords)[gmx_index_phase[phase][atom]][(itim->normal+1)%3];
			   y = (gmx_coords)[gmx_index_phase[phase][atom]][(itim->normal+2)%3];
			   z = (gmx_coords)[gmx_index_phase[phase][atom]][itim->normal];
			   if(phase==SUPPORT_PHASE)radius = itim->radii[atom];

                 	   itim->phase[phase][3*atom] = x;
                 	   itim->phase[phase][3*atom+1] = y;
                 	   itim->phase[phase][3*atom+2] = z;
                 	   itim->phase_index[phase][atom] =  atom; 
		           for(sign=-1; sign<=1; sign+=2){ 
			       if(sign*x > itim->box[0]/2.- scale*itim->alpha) {  /*note that -x > box - a <=>  x < - (box-a)  */
			            itim->phase[phase][3*(nindex[phase]+additional)] = x-sign*itim->box[0];
			            itim->phase[phase][3*(nindex[phase]+additional)+1] = y;
			            itim->phase[phase][3*(nindex[phase]+additional)+2] = z;
                 	   	    itim->phase_index[phase][nindex[phase]+additional] =  atom; 
                 	   	    if(phase==SUPPORT_PHASE)itim->radii[nindex[phase]+additional] =  radius; 
			            additional++;
			       }	
			       if(sign*y > itim->box[1]/2.- scale*itim->alpha) {
			            itim->phase[phase][3*(nindex[phase]+additional)] = x;
			            itim->phase[phase][3*(nindex[phase]+additional)+1] = y-sign*itim->box[1];
			            itim->phase[phase][3*(nindex[phase]+additional)+2] = z;
                 	   	    itim->phase_index[phase][nindex[phase]+additional] =  atom; 
                 	   	    if(phase==SUPPORT_PHASE)itim->radii[nindex[phase]+additional] =  radius; 
			            additional++;
		                    for(sign2=-1; sign2<=1; sign2+=2){ 
			                 if(sign2*x > itim->box[0]/2.- scale*itim->alpha) {
			                      itim->phase[phase][3*(nindex[phase]+additional)] = x-sign2*itim->box[0];
			                      itim->phase[phase][3*(nindex[phase]+additional)+1] = y-sign*itim->box[1];
			                      itim->phase[phase][3*(nindex[phase]+additional)+2] = z;
                 	   	    	      itim->phase_index[phase][nindex[phase]+additional] =  atom; 
                 	   	    	      if(phase==SUPPORT_PHASE)itim->radii[nindex[phase]+additional] =  radius; 
			                      additional++;
			                 }	
			            }
			       }	
			    }
			    /* if the memory allocated so far is not enough, let's get another bunch...*/
			    if(additional + nindex[phase] > realloc_factor*nindex[phase]-8) { /* 8 because next round we can add at most 8 new periodic copies of one atom */
					realloc_factor*=2;
                                        itim->phase[phase] = (real * ) realloc(itim->phase[phase],realloc_factor*nindex[phase]*sizeof(real) * 3);
                                        itim->phase_index[phase] = (int * ) realloc(itim->phase_index[phase],realloc_factor*nindex[phase]*sizeof(int) );
            		                if(phase==SUPPORT_PHASE)itim->radii=  (real *)realloc(itim->radii, realloc_factor*nindex[phase]*sizeof(real));
			    }

			}
		        itim->n[phase] = nindex[phase] + additional;	
            		if(phase==SUPPORT_PHASE)itim->radii=  (real *)realloc(itim->radii,itim->n[phase]*sizeof(real));
                        itim->phase[phase] = (real * ) realloc(itim->phase[phase],itim->n[phase]*sizeof(real) * 3);
                        itim->phase_index[phase] = (int * ) realloc(itim->phase_index[phase],itim->n[phase]*sizeof(int) );
			break;
                        default: exit(printf("Geometry not implemented [arrange_datapoints()]\n"));
		  }	
	
		break;
	        default: exit(printf("Periodicity not implemented\n")); break;
	  } // end switch periodicity
	}
        if(itim->nphases>itim->ngmxphases) itim->n[itim->nphases-1] = itim->n[INNER_PHASE]+itim->n[OUTER_PHASE];
        itim->mask=mask;
}

Histogram * histo_init(Histogram * histo, int N, int nbins, double range) { 

	histo = (Histogram*)malloc(N*sizeof(Histogram));
	if(histo==NULL) exit(printf("Internal error\n"));
	histo->N=N;
	histo->nbins= nbins; // we are sampling half of the boxlength
	histo->iterations= 0;
	histo->rdata = (double *)calloc(histo->N*histo->nbins,sizeof(double));
	histo->size = histo->minsize = range;
	histo->bWidth = histo->size / histo->nbins;
	return  histo;
}

static int populate_histogram(int index, real dist, Histogram * histo, ITIM * itim, real value){
            /* we have values from -box/2 to box/2, let's shift them back to gromacs convention (0:box)*/
	    dist += histo->size/2.;
	    if(dist<histo->size && dist>=0){ 
	 		return add_histo(histo,index,(double)dist,(double)value);
	    } 	
	    return 0;
}

void check_inclusive_conditions(atom_id ** index, int *gnx){

	int i,c=0;
        const char msg[]="\nConditions for the -inclusive options are not fulfilled\n";		
	if ( index[INNER_PHASE][0] < index[OUTER_PHASE][0] )   { 
		for(i=0;i<gnx[INNER_PHASE];i++,c++)
			if(index[INNER_PHASE][i] != index[SUPPORT_PHASE][c] ) exit(printf(msg));
		for(i=0;i<gnx[OUTER_PHASE];i++,c++)
			if(index[OUTER_PHASE][i] != index[SUPPORT_PHASE][c] ) exit(printf(msg));
	} else { 
		for(i=0;i<gnx[OUTER_PHASE];i++,c++)
			if(index[OUTER_PHASE][i] != index[SUPPORT_PHASE][c] ) exit(printf(msg));
		for(i=0;i<gnx[INNER_PHASE];i++,c++)
			if(index[INNER_PHASE][i] != index[SUPPORT_PHASE][c] ) exit(printf(msg));
	}

}



ITIM * init_intrinsic_surface(Flags myflags, int normal, real alpha, real mesh,  matrix box, int ngrps, int *nbins, int maxlayers, real * radii, int** gmx_index,int *gnx, int *com_opt,int dump_mol, const char ** geometry, int ngrps_add, t_topology * top){
	    ITIM * itim;
	    int i;
            int histofactor=1;
            /* ngrps + 1 here  because of the random phase used to to a MC estimate of the bin volumes for the density profile*/
     	    init_itim(ngrps+1,ngrps_add) ;
	    itim = global_itim;
	    itim->info = myflags.bInfo; // TODO put me into cmdline...
	    itim->nadd_index =  ngrps_add;
	    itim->maxlayers = maxlayers;
     	    itim->method = METHOD_ITIM;
            itim->bMol=myflags.bMol;
            itim->bOrder=myflags.bOrder;
            itim->bVirial=myflags.bVirial;
            itim->bInclusive=myflags.bInclusive;
            itim->bMCnormalization=myflags.bMCnormalization;
            itim->n_histo=GET_HISTO_SIZE();
            itim->dump_mol=dump_mol;
            /* 1 mass dens x phases  (+random phase),  ( 1 number dens , 4 order, 4 error on order) x phases */
 //TODO: comment this better
            if(itim->bOrder) exit(printf("Internal error: this has to be redone\n"));
     	    itim->alpha=alpha;
	    itim->skin=0.0;  // TODO: CMDLINE ?
	    itim->range=10.0; // TODO: CMDLINE ? check the other comment about range...
	    switch(geometry[0][0]){
	    	case 'p':  itim->geometry = SURFACE_PLANE; 
     	                   itim->normal=normal; 
      	                   if(itim->method==METHOD_ITIM)  itim->periodic[SUPPORT_PHASE] = NONE; 
			   break;
		default:  exit(printf("Geometry not implemented so far.\n"));
	    }

	    itim->target_mesh_size= mesh;
            /* itim->com_opt[] are by default zero */
            for(i=0;i<itim->ngmxphases;i++)
	         itim->com_opt[i] = com_opt[i];
	    itim->histograms = histo_init(itim->histograms, itim->n_histo, *nbins,(double) box[itim->normal][itim->normal]) ; 
#ifdef INTEGER_HISTO
            int_histo = malloc(itim->n_histo*(*nbins)*sizeof(unsigned long long int));
            for(i=0;i<itim->n_histo*(*nbins);i++) int_histo[i]=0; 
#endif
            itim->radii=  (real *)malloc(gnx[SUPPORT_PHASE]*sizeof(real));
	    for(i=0;i<gnx[SUPPORT_PHASE];i++){
		itim->radii[i] = radii[gmx_index[SUPPORT_PHASE][i]];
	    }
            itim->gmx_index = gmx_index;
	    if(itim->bMol) {  // TODO: check if some of the other spol_atom2molindex() calls can be replaced by referencing this index here...
		int isizem=gnx[SUPPORT_PHASE];  
		itim->indexm = (int*) malloc(isizem*sizeof(int));
		itim->backindex= (int*) malloc((top->mols.nr+1)*sizeof(int));
		for(i=0;i<gnx[SUPPORT_PHASE];i++) itim->indexm[i]=itim->gmx_index[SUPPORT_PHASE][i]; // NOTE: TODO check "additional"  under the PATCH case
                spol_atom2molindex(&isizem, itim->indexm,itim->backindex, &(top->mols));
            }
	    return itim;
}

void finalize_intrinsic_profile(real *** density, int * nslices, real * slWidth){
	int i,j;
	ITIM *itim = global_itim;
	Histogram * histo = itim->histograms;
  	*density = (real **) malloc( itim->maxlayers*itim->n_histo* sizeof(real*));
	for(i=0;i<itim->n_histo;i++){
	    (*density)[i]=(real*) malloc(histo->nbins * sizeof(real));
        }
	/* Let's re-determine the number of slices, given the minimum box-size found during the run*/
	*nslices = (int)((histo->minsize*(histo->nbins-1)) / histo->size);  
	*slWidth = histo->bWidth ;
    	for(i=0;i<itim->n_histo;i++){
  	  for(j=0;j<histo->nbins;j++){
#ifndef TO_INTEGRATE_CHARGES
  		histo->rdata[ i * histo->nbins + j ] /= (histo->iterations * *slWidth); 
#endif
                switch(itim->geometry){
  		  case SURFACE_PLANE: 
#ifndef TO_INTEGRATE_CHARGES
                             histo->rdata[ i * histo->nbins + j ]/=(2*(itim->box[0] * itim->box[1]));
			     //printf("CALLED i=%d j=%d h=%f\n",i,j, histo->rdata[ i * histo->nbins + j ]);
#endif
  			  break;
                  case SURFACE_SPHERE: 
                  case SURFACE_GENERIC: 
  			break;
                  default: break;
                 }
  	   }
	}
        /* let's apply normalizations...*/
        if(itim->bMCnormalization){
    	   for(j=SUPPORT_PHASE ; j<itim->nphases+itim->nadd_index;j++){
	      if (j==itim->RANDOM_PHASE) continue;
	      for(i=0;i<itim->maxlayers+1;i++){
	         if(j==SUPPORT_PHASE && i==0) continue; 
	         if(j!=SUPPORT_PHASE && i!=0) continue;
	         for(int m=ATOMIC;m<=MOLECULAR;m++){
      		    if (itim->bMol==FALSE && m==MOLECULAR) continue;
		    if(j!=SUPPORT_PHASE && m==MOLECULAR) continue;
	            int index     = GET_HISTO_INDEX(INTRINSIC_DENSITY,j,                 i,m     ,__LINE__);
	            int index_rnd = GET_HISTO_INDEX(INTRINSIC_DENSITY,itim->RANDOM_PHASE,0,ATOMIC,__LINE__);
		    if(index==-1000 || index_rnd == -1000) continue;
		    printf("Normalizing index %d\n",index);
  	            for(int k=0;k<histo->nbins;k++){
                         double norm =  histo->rdata[index_rnd*histo->nbins + k ];
			// printf("Normalizing index=%d bin %d using norm %g\n",index,k,norm);
                         if(norm>0) histo->rdata[ index * histo->nbins + k ] /= norm ;
                    }
                 }
              }
           }
        }
  
  	for(i=0;i<itim->n_histo;i++){
  	   for(j=0;j<histo->nbins;j++){
  		(*density)[i][j] = (real)histo->rdata[ i * histo->nbins+j];
             }
  	}
}
int compute_intrinsic_surface(int bCluster, matrix box, int ngrps, rvec * gmx_coords, int *nindex, atom_id ** gmx_index_phase,t_topology * top, int natoms, t_pbc pbc){
	
	ITIM * itim = global_itim;
	int ind=itim->normal;
        static int * mask=NULL;
	int success=1;
#ifdef TIME_PROFILE
        struct timeval tp;
        struct timeval tp2;
#endif

	ind=(ind+1)%3;  itim->box[0] = box[ind][ind];
	ind=(ind+1)%3;  itim->box[1] = box[ind][ind];
	ind=(ind+1)%3;  itim->box[2] = box[ind][ind];
#ifdef TIME_PROFILE
        gettimeofday(&tp, NULL);
#endif
	switch (itim->method) { 
		
		case METHOD_ITIM: 
	              /* creates a grid which is as close as possible to the target mesh size, 
	                 but still has an integer number of cells in the sim box*/
	              init_itim_grid(itim);
		      /* rearrange positions and indices to be fed to itim routines */
		      /* New mask: -1 -> not in the main cluster
				    0 -> layer not assigned
				    i -> i-th layer (positive or negative) */ 
                      generate_mask_ns(bCluster,itim->bInclusive,gmx_coords,&mask,nindex, (int**)gmx_index_phase,box,itim->ngmxphases,natoms,pbc); 
	              arrange_datapoints(itim, gmx_coords, nindex,  (int**)gmx_index_phase,mask);
                      reinit_mask(itim);
		      for(itim->current_layer = 1 ; itim->current_layer <= itim->maxlayers  ; itim->current_layer ++){
		      	success *= compute_itim_points(DIR_POSITIVE,itim,(int**)gmx_index_phase,top,gmx_coords,mask);
		      }
		      kd_free(itim->mesh.tree); itim->mesh.tree=NULL;
		break;
		default: exit(printf("Method not implemented yet\n"));
		break;
	}
#ifdef TIME_PROFILE
        gettimeofday(&tp2, NULL);
        fprintf(stdout,"Time to build surface (method: %s): millisec=%f\n",itim->method_name[itim->method],1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
	if (success==0) {
	   // reset all lines and return 0, the calling function should decide what to do in this case
           // TODO BUG: this still does not fix the problem, after one "wrong" frame identification keeps failing
	   for(int i=0;i<itim->mesh.nelem;i++){
             	itim->mesh.flag[i]=DIR_POSITIVE;   
           }
	   return 0;
        }
        arrange_alpha_points (itim,(int**)gmx_index_phase,top,gmx_coords); 
	if(itim->info)fprintf(stdout,"Number of surface elements = %d\n",itim->nalphapoints);
	return 1;
}

void compute_histogram(matrix box,atom_id ** gmx_index_phase,t_topology * top, char dens_opt, t_trxframe * fr){  
/*************************************

          TODO TODO TODO 

This is too cluttered. Reorganize the code...

          TODO TODO TODO 

 ************************************/
    ITIM * itim=global_itim;
    Histogram *histo = itim->histograms; 
    static real * result_points = NULL;
    real dist =0,check=0.0;
    int i,j,k,dim,sampled;
    enum { QSORT , KDTREE} SORT_METHOD  ;
    real *ppoints=NULL;
    real size,size2;
    real order[4]={0.,0.,0.,0.};
    real order2[4]={0.,0.,0.,0.};
    real minbox;
    int pn1=0;
    rvec normal;
    SORT_METHOD = KDTREE;
    real p4[3];
    static int rp_n=-1;
    real tmpreal;
    rvec tmprvec;
#ifdef TIME_PROFILE
    struct timeval tp;
    struct timeval tp2;
    gettimeofday(&tp, NULL);
#endif
    if(result_points==NULL) { 
		rp_n = 3*itim->nalphapoints;
		result_points = (real *) realloc(result_points,rp_n * sizeof(real)); // TODO: avoid creating this array, just use the tree
    }
    if(rp_n==-1) exit(printf("Internal error, rp_n not assigned correctly in compute_histogram()\n"));
    switch(itim->geometry){
	case SURFACE_PLANE:   size=itim->box[2]/2; size2=size*size; break;
        case SURFACE_SPHERE: 
        case SURFACE_GENERIC: 
			      size = my_min(my_min(itim->box[2],itim->box[1]),itim->box[0])/2; size2=size*size; break;
        default: exit(printf("Not implemented\n"));
    }

    if(2*size<histo->minsize) histo->minsize=2*size;
    histo->iterations++;
    for(j=SUPPORT_PHASE ; j<itim->nphases;j++){

        if(j==itim->RANDOM_PHASE && !itim->bMCnormalization) continue;

	int isizem, *indexm,*backindex;
	int molecular_layer,atomic_layer;
        if(itim->bMol) { 
		isizem=itim->n[j];
		snew(indexm,isizem);
		snew(backindex,top->mols.nr+1);
	}


	if(!(j==itim->RANDOM_PHASE && itim->bMCnormalization) && itim->bMol ) {
		for(i=0;i<itim->n[j];i++) indexm[i]=itim->gmx_index[j][i]; // NOTE: TODO check "additional"  under the PATCH case
                spol_atom2molindex(&isizem, indexm,backindex, &(top->mols));
	}
  	/*This is Miguel's original algorithm. */
	for(i=0;i<itim->n[j];i++){
		molecular_layer=0;
		atomic_layer=0;
                real locmass=0.0,locm,locnumber;
                rvec v1,v2,vm,vn,newpos,oldpos;
		if(j==SUPPORT_PHASE){
			if(itim->bMol) { 
			   int mol = indexm[i];
			   int i1 = backindex[mol] ;
			   int i2 = backindex[mol+1];
			   if(i2==-1) i2=itim->n[SUPPORT_PHASE]; // i2==-1 in case of last atom.
			   int minlayer=123456;
			   int k;
			   for(k=i1;k<i2;k++){
			   	if (itim->mask[k] >0 && itim->mask[k] < minlayer ) { 
			   		minlayer = itim->mask[k];
			   	}
			   }	
			   // with the molecular layer we add all atoms belonging to the same molecule into 
                           // the outermost layer associated to any of the molecule's atoms. Clear, uh?
			   molecular_layer= ( (minlayer==123456) ? 0 :minlayer) ;	
			}
			atomic_layer =itim->mask[i]; 
		}
                if(itim->com_opt[j]) {
                     if(j>=itim->ngmxphases){  exit(printf("Internal error\n")) ;} 
                     p4[0]=p4[1]=p4[2]=0.0;
                     for(k=0;k<itim->com_opt[j];k++){ // let's go through the molecule's atoms
                            locm = itim->masses[gmx_index_phase[j][i+k]];
			    locmass += locm;
                            for(dim=0;dim<3;dim++){ 
                                newpos[dim]=itim->phase[j][3*(i+k)+dim] ;
                                if(k>0){// remove pbc...
                                     while(newpos[dim]-oldpos[dim] >  itim->box[dim]/2.) {newpos[dim]-= itim->box[dim]; }
                                     while(newpos[dim]-oldpos[dim] < -itim->box[dim]/2.) {newpos[dim]+= itim->box[dim]; }
                                }
                                p4[dim]+= locm * newpos[dim] ; // and assign the point in the phase...
                                oldpos[dim]=newpos[dim];
                            }
                     }
                     if (itim->bOrder){ 
                                for(dim=0;dim<3;dim++){
                                   vm[dim]=itim->phase[j][3*(i+1)+dim];
                                   vn[dim]=itim->phase[j][3*(i)+dim];
                                }
                                rvec_sub(vm,vn,v1); // H1-O 
                                for(dim=0;dim<3;dim++){
                                   vm[dim]=itim->phase[j][3*(i+2)+dim];
                                }
                                rvec_sub(vm,vn,v2); // H1-O 
                                rvec_add(v1,v2,vm);
                                unitv(vm,vm);
                                cprod(v1,v2,vn);
                                unitv(vn,vn);
                     }
                     i+=itim->com_opt[j]-1;
                     p4[0]/=locmass; p4[1]/=locmass; p4[2]/=locmass;
                } else { 
                     if(itim->bOrder) { order[0]=order[1]=order[2]=order[3]=0.0; 
                     		        order2[0]=order2[1]=order2[2]=order2[3]=0.0; } 
                     if(j<itim->ngmxphases){ /* to take in account the random phase */
                            real value;
		            real vec[3],wec[3],Vec[3],tmpv,mod;
                            switch(dens_opt){
			    	case 'm': value = itim->masses[gmx_index_phase[j][i]]; break;
			    	case 'c': value = itim->charges[gmx_index_phase[j][i]]; break;
			    	case 'n': value = 1; break;
			    	case 'x': value = 0.5 *itim->masses[gmx_index_phase[j][i]]* SQR(fr->v[gmx_index_phase[j][i]][0]); break;
			    	case 'y': value = 0.5 *itim->masses[gmx_index_phase[j][i]]* SQR(fr->v[gmx_index_phase[j][i]][1]); break;
			    	case 'z': value = 0.5 *itim->masses[gmx_index_phase[j][i]]* SQR(fr->v[gmx_index_phase[j][i]][2]); break;
			    	case 'X': value = fr->v[gmx_index_phase[j][i]][0]; break;
			    	case 'Y': value = fr->v[gmx_index_phase[j][i]][1]; break;
			    	case 'Z': value = fr->v[gmx_index_phase[j][i]][2]; break;
			    	case 'H': if((i%3)==0) { value=0; break; } 
					  mod=0;
					  for(int c=0;c<3;c++){
					     vec[c]=0.0;
					     tmpv=fr->x[gmx_index_phase[j][3*(i/3)+1]][c]-
					            fr->x[gmx_index_phase[j][3*(i/3)]][c];
				             while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				             while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					     vec[c]+=tmpv;

					     tmpv=fr->x[gmx_index_phase[j][3*(i/3)+2]][c]-
					            fr->x[gmx_index_phase[j][3*(i/3)]][c];
				             while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				             while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					     vec[c]+=tmpv;

					     mod+=vec[c]*vec[c];
				   	  }
				  	  mod=sqrt(mod);
					  for(int c=0;c<3;c++) vec[c]/=mod;
					  value = 3*vec[2]*vec[2] - 1  ;
					  break;

				case 'h': if((i%3)==0) { value=0; break; } 
					  mod=0;
					  for(int c=0;c<3;c++){
					    Vec[c]=0.0;
					    tmpv=fr->x[gmx_index_phase[j][3*(i/3)+1]][c]-
					           fr->x[gmx_index_phase[j][3*(i/3)]][c];
				            while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				            while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					    vec[c]=tmpv;

					    tmpv=fr->x[gmx_index_phase[j][3*(i/3)+2]][c]-
					           fr->x[gmx_index_phase[j][3*(i/3)]][c];
				            while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				            while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					    wec[c]=tmpv;
					  }
                                	  cprod(vec,wec,Vec);
                                          unitv(Vec,Vec);
					  value = 3*Vec[2]*Vec[2] - 1  ;
					  break;

				case 'o': if((i%3)!=0) { value=0; break; } 
					  mod=0;
					  for(int c=0;c<3;c++){
					    Vec[c]=0.0;
					    tmpv=fr->x[gmx_index_phase[j][3*(i/3)+1]][c]-
					           fr->x[gmx_index_phase[j][3*(i/3)]][c];
				            while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				            while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					    vec[c]=tmpv;

					    tmpv=fr->x[gmx_index_phase[j][3*(i/3)+2]][c]-
					           fr->x[gmx_index_phase[j][3*(i/3)]][c];
				            while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				            while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					    wec[c]=tmpv;
					  }
                                	  cprod(vec,wec,Vec);
                                          unitv(Vec,Vec);
					  value = 3*Vec[2]*Vec[2] - 1  ;
					  break;



				case 'O': if((i%3)!=0) {value=0; break;}
					  mod=0;
					  for(int c=0;c<3;c++){
					     vec[c]=0.0;
					     tmpv=fr->x[gmx_index_phase[j][3*(i/3)+1]][c]-
					            fr->x[gmx_index_phase[j][3*(i/3)]][c];
				             while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				             while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					     vec[c]+=tmpv;

					     tmpv=fr->x[gmx_index_phase[j][3*(i/3)+2]][c]-
					            fr->x[gmx_index_phase[j][3*(i/3)]][c];
				             while(tmpv>itim->box[c]/2.) tmpv-=itim->box[c];
				             while(tmpv<-itim->box[c]/2.) tmpv+=itim->box[c];
 					     vec[c]+=tmpv;

					     mod+=vec[c]*vec[c];
				   	  }
				  	  mod=sqrt(mod);
					  for(int c=0;c<3;c++) vec[c]/=mod;
					  value = 3*vec[2]*vec[2]-1;
					  break;
				case 'K': if((i%3)!=0) {value=0; break;}
					  mod=0;
					  for(int c=0;c<3;c++){
					    mod+= 0.5*itim->masses[gmx_index_phase[j][i+0]] * SQR(fr->v[gmx_index_phase[j][i+0]][c]);
					    mod+= 0.5*itim->masses[gmx_index_phase[j][i+1]] * SQR(fr->v[gmx_index_phase[j][i+1]][c]);
					    mod+= 0.5*itim->masses[gmx_index_phase[j][i+2]] * SQR(fr->v[gmx_index_phase[j][i+2]][c]);
				   	  }
					  value = mod;
					  break;

			
#ifdef VIRIAL_EXTENSION
			    	case 't': 
					tmpreal = itim->box[0]*itim->box[1]*itim->box[2];

					if(itim->bVirial==true){
						copy_rvec(fr->vir[gmx_index_phase[j][i]],tmprvec);
						tmpreal = itim->box[0]*itim->box[1]*itim->box[2];
						value = (0.5 * (tmprvec[0]+tmprvec[1]));  
					} else { 
                                                real mass = itim->masses[gmx_index_phase[j][i]] ;

						/* 
						   == For testing purposes ==
						note: 1 kJ/mol = 16.6054 bar nm^3 

						The following should give the total kinetic energy in kJ/mol :
					  	   copy_rvec(fr->v[gmx_index_phase[j][i]],tmprvec);
						   value =  0.5 * mass * (SQR(tmprvec[2]) + SQR(tmprvec[1])+SQR(tmprvec[0])) ;

                                                The following should give the total vir_xx in bar nm^3 (mind the sign & factor)
						   copy_rvec(fr->vir[gmx_index_phase[j][i]],tmprvec);
                                                   value = -tmprvec[0]/2.  ; 
						*/

						/* we write out the average (Pxx+Pyy)/2. */
						// ideal gas 
						copy_rvec(fr->v[gmx_index_phase[j][i]],tmprvec);
						value =  mass * 0.5 * (SQR(tmprvec[0])+ SQR(tmprvec[1])) * 16.6054 ; 
						// virial (the factor -2 is already in)
						copy_rvec(fr->vir[gmx_index_phase[j][i]],tmprvec);
                                                value += 0.5 * (tmprvec[0] + tmprvec[1])  ; 
					}
						break ; 
			    	case 'E': 
						tmpreal = fr->pener[gmx_index_phase[j][i]];
						value = tmpreal / (itim->box[0]*itim->box[1]*itim->box[2]);
						break ; 
			    	case 'U': 
						tmpreal = 0;
						for(int ii=0;ii<itim->n[j];ii++)
                                                       tmpreal+=fr->pener[gmx_index_phase[j][ii]];
						value = tmpreal / (itim->box[0]*itim->box[1]*itim->box[2]);
						break ; 
#endif
			    	default : value =1 ; break;
                            }
                            locmass=value;
		            p4[0]=itim->phase[j][3*i]; p4[1]=itim->phase[j][3*i+1]; p4[2]=itim->phase[j][3*i+2];
                     } else { 
                            locmass=itim->box[0]*itim->box[1]*itim->box[2]/itim->n[j];
			    p4[0]=(((real)rand()/RAND_MAX)-0.5)*itim->box[0];
			    p4[1]=(((real)rand()/RAND_MAX)-0.5)*itim->box[1];
			    p4[2]=(((real)rand()/RAND_MAX)-0.5)*itim->box[2];
                     }
                }
#ifdef GITIM_LR
		if(NULL!=getenv("GITIM_LEFT")  && p4[2]>0 ) continue;
		if(NULL!=getenv("GITIM_RIGHT")  && p4[2]<0 ) continue;
#endif
                switch(itim->geometry){
                                case SURFACE_PLANE: if(fabs(p4[2])>size) continue; break;
                                default: exit(printf("No such a geometry (%d)\n",itim->geometry));
                }

                if(itim->geometry==SURFACE_PLANE && j==INNER_PHASE)
	             if(fabs(p4[0])>=itim->box[0]/2. ||  fabs(p4[1]) >=  itim->box[1]/2. || fabs(p4[2])>=itim->box[2]/2.) 
                           continue; 
		 global_real_pointer = p4;
                 if(itim->geometry==SURFACE_PLANE && SORT_METHOD==QSORT) {
     		        qsort((void*)itim->alphapoints, itim->nalphapoints, 3*sizeof(real), periodic_distance );
		        ppoints=itim->alphapoints; 
        	        pn1=itim->nalphapoints ; 
                 } else { // Not using Qsort but kdtrees
                        switch(itim->geometry){
			   case SURFACE_PLANE:
/* TODO: implement local normal */
		 	       dist = perform_interpolation(Surface, p4, &(normal[0]),itim,i,j);
			   break;
                           default: exit(printf("No such a geometry (%d)\n",itim->geometry));
                        }
		 }
		 if(itim->bOrder && itim->com_opt[j]){
			        exit(printf("Re-implement the calculation of the normal vector in interpolate_distance3D_2() \n"));
                                switch(itim->geometry){
					case SURFACE_PLANE:
                                     		v1[0]=0;v1[1]=0;v1[2]=1.0;
                                     		v2[0]=normal[0];v2[1]=normal[1];v2[2]=normal[2];
                                              break;
                                        default: exit(printf("Surface %d not implemented\n",itim->geometry)) ; 
                                              break;
                                }
				order[0]=iprod(v1,vm);                                         order2[0]=order[0]*order[0];
				order[1]=iprod(v1,vn); order[1]=(3*order[1]*order[1] -1. )/2.; order2[1]=order[1]*order[1];
				order[2]=iprod(v2,vm);                                         order2[2]=order[2]*order[2];
				order[3]=iprod(v2,vn); order[3]=(3*order[3]*order[3] -1. )/2.; order2[3]=order[3]*order[3];
                 }
                 locnumber=1;
#ifdef TO_INTEGRATE_CHARGES
		 sampled=populate_histogram(j, dist, histo, itim,(real)((int)(1000*locmass)));
#else
		 if(j==SUPPORT_PHASE){
			int add,phase_index;
			for(add=-1;add<itim->nadd_index;add++){
			   if(add==-1){ 
				phase_index=SUPPORT_PHASE; 
			   } else {
				phase_index=itim->ngmxphases+add;	
				if(itim->mask_add[add][i]==0) continue;
			   }
			   if(atomic_layer>0){
				//printf("phase_index = %d\n",phase_index);
		 	   	sampled=populate_histogram(GET_HISTO_INDEX(INTRINSIC_DENSITY,phase_index,atomic_layer,ATOMIC,__LINE__), dist, histo, itim,(real)(locmass));
		 	   	sampled=populate_histogram(GET_HISTO_INDEX(LAYER_DISTRIBUTION,phase_index,atomic_layer,ATOMIC,__LINE__), p4[2], histo, itim,2.*(real)(locmass));
			   }
			   if(molecular_layer>0 && itim->bMol){
			   	sampled=populate_histogram(GET_HISTO_INDEX(INTRINSIC_DENSITY,phase_index,molecular_layer,MOLECULAR,__LINE__), dist, histo, itim,(real)(locmass));
			   	sampled=populate_histogram(GET_HISTO_INDEX(LAYER_DISTRIBUTION,phase_index,molecular_layer,MOLECULAR,__LINE__), p4[2], histo, itim,2*(real)(locmass));
			   }
			}
		 } else { 
		 	   sampled=populate_histogram(GET_HISTO_INDEX(INTRINSIC_DENSITY,j,0,ATOMIC,__LINE__), dist, histo, itim,(real)(locmass));
		 	   sampled=populate_histogram(GET_HISTO_INDEX(LAYER_DISTRIBUTION,j,0,ATOMIC,__LINE__),p4[2], histo, itim,2*(real)(locmass));
		 }
#endif
                 //printf("SAMPLING %d %d %f %d\n",i,sampled,locmass,check+=locmass);
                 //printf("SAMPLING %d %d\n",i,sampled);
               
#if 0
                 if(normal[0]!=NORMAL_UNDEFINED) /*SAW: TODO NOTE that this creates an inconsistency between the counting of masss profile and the others. Should one just drop points which do not have a triangle associated on the surface, or should we use the macroscopic vector in those cases? See perform_interpolation() for the handling of pathological cases.*/
                    if( j < itim->ngmxphases ) { 
		              populate_histogram(OFF_NUMBER+j, dist, histo, itim,locnumber);
                              if(itim->bOrder && normal[0]!=NORMAL_UNDEFINED ) { 
			         populate_histogram(OFF_ORDER1+j, dist, histo, itim,  order[0]);
			         populate_histogram(OFF_ORDER1_2+j, dist, histo, itim,order2[0]);
			         populate_histogram(OFF_ORDER2+j, dist, histo, itim,  order[1]);
				 populate_histogram(OFF_ORDER2_2+j, dist, histo, itim,order2[1]);
			         populate_histogram(OFF_ORDER3+j, dist, histo, itim,  order[2]);
			         populate_histogram(OFF_ORDER3_2+j, dist, histo, itim,order2[2]);
			         populate_histogram(OFF_ORDER4+j, dist, histo, itim,  order[3]);
			         populate_histogram(OFF_ORDER4_2+j, dist, histo, itim,order2[3]);
                             }
                    }
#endif
	}
	if(itim->bMol) { 
	   free(indexm);
	   free(backindex);
	}
    }
#ifdef TIME_PROFILE
    gettimeofday(&tp2, NULL);
    fprintf(stdout,"Time to make the histogram(s) for %d phases: millisec=%f\n",itim->nphases,1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
}

void reset_counters(){
	kd_free(Surface);
	Surface=NULL;
	ITIM * itim = global_itim;
	itim->mesh.nelem=0;
	itim->nalphapoints = 0;
}

void  compute_intrinsic_profile(matrix box, atom_id ** gmx_index_phase, t_topology * top, char dens_opt, t_trxframe * fr){
	compute_histogram(box,gmx_index_phase,top,dens_opt,fr);
	reset_counters();
}



/**/

typedef struct {
  char *atomname;
  int nr_el;
} t_electron;

/****************************************************************************/
/* This program calculates the partial density across the box.              */
/* Peter Tieleman, Mei 1995                                                 */
/****************************************************************************/

/* used for sorting the list */
int compare(void *a, void *b)
{
  t_electron *tmp1,*tmp2;
  tmp1 = (t_electron *)a; tmp2 = (t_electron *)b;

  return strcmp(tmp1->atomname,tmp2->atomname);
}

int get_electrons(t_electron **eltab, const char *fn)
{
  char buffer[256];  /* to read in a line   */
  char tempname[80]; /* buffer to hold name */
  int tempnr; 

  FILE *in;
  int nr;            /* number of atomstypes to read */
  int i;

  if ( !(in = gmx_ffopen(fn,"r")))
    gmx_fatal(FARGS,"Couldn't open %s. Exiting.\n",fn);

  if(NULL==fgets(buffer, 255, in))
  {
      gmx_fatal(FARGS,"Error reading from file %s",fn);
  }
 
  if (sscanf(buffer, "%d", &nr) != 1)
    gmx_fatal(FARGS,"Invalid number of atomtypes in datafile\n");

  snew(*eltab,nr);

  for (i=0;i<nr;i++) {
    if (fgets(buffer, 255, in) == NULL)
      gmx_fatal(FARGS,"reading datafile. Check your datafile.\n");
    if (sscanf(buffer, "%s = %d", tempname, &tempnr) != 2)
      gmx_fatal(FARGS,"Invalid line in datafile at line %d\n",i+1);
    (*eltab)[i].nr_el = tempnr;
    sprintf((*eltab)[i].atomname ,"%s",tempname);
  }
  gmx_ffclose(in);
  
  /* sort the list */
  fprintf(stdout,"Sorting list..\n");
  qsort ((void*)*eltab, nr, sizeof(t_electron), 
	 (int(*)(const void*, const void*))compare);

  return nr;
}

void remove_phase_pbc(int ePBC, t_atoms *atoms, matrix box, rvec x0[], int axis, atom_id *index, int index_nr){
     /*we assume here that atoms have been already put into the box */
     /* compute the density at different control points: box edges, middle + some more */
     int rho[5],rho_max,i,nbins=25; 
     static real shift=0.0;
     real bWidth ,z ;
     ITIM* itim = global_itim;
  // TODO: check: this does not work with a solid ... bin must bigger than average interparticle z distance...
     bWidth = box[axis][axis]/nbins;
     while (1){ 
	 if(shift!=0.0){
     	   for(i=0; (i<atoms->nr); i++) {
	   	x0[i][axis]+=shift;
	   }
	   put_atoms_in_box(ePBC,box,atoms->nr,x0);
         }
         rho[0]=rho[1]=rho[2]=rho[3]=rho[4]=0;
         for(i=0; (i<index_nr); i++) {
             z=x0[index[i]][axis];

             if(z<bWidth) rho[0]++; 
             else if (z>box[axis][axis]-bWidth) rho[4]++;
             else if (z>(nbins/2)*bWidth && z < (nbins/2+1)*bWidth) rho[2]++;
             else if (z>(nbins/4)*bWidth &&  z < (nbins/4+1)*bWidth) rho[1]++;
             else if (z>(3*nbins/4)*bWidth && z < (3*nbins/4+1)*bWidth) rho[3]++;

         }
         for(rho_max=rho[0],i=1; i<5;i++) if (rho[i]>rho_max) rho_max=rho[i];
         if(rho[0] > rho_max/2 || rho[4] > rho_max/2) { /* careful: not >= otherwise the algorithm doesn't work when rho=0 
            					       in all sampled regions */
            	shift+=bWidth*rand()/RAND_MAX;
		if(itim->info) fprintf(stdout,"Trying to shift the box by %f nm\n",shift);
		if(shift>box[axis][axis]) exit(printf("Error: it was not possible to center the phase in the box\n"));
         } else return ;
    }
}

void center_coords(t_atoms *atoms,matrix box,rvec x0[],int axis,atom_id *index, int index_nr)
{
  int  i,m;
  real tmass,mm;
  rvec com,shift,box_center;
  
  tmass = 0;
  clear_rvec(com);
  clear_rvec(shift);

  if(index==NULL) { 
            for(i=0; (i<atoms->nr); i++) {
              mm     = atoms->atom[i].m;
              tmass += mm;
              for(m=0; (m<DIM); m++) 
                com[m] += mm*x0[i][m];
            }
  } else { 
            for(i=0; (i<index_nr); i++) {
              mm     = atoms->atom[index[i]].m;
              tmass += mm;
              for(m=0; (m<DIM); m++) 
                com[m] += mm*x0[index[i]][m];
            }
  }
  for(m=0; (m<DIM); m++) 
    com[m] /= tmass;
  calc_box_center(ecenterDEF,box,box_center);
  rvec_sub(box_center,com,shift);
  
  if(index==NULL){
      shift[axis] = box_center[axis];
  } else { 
    for(m=0; (m<DIM) ; m++) shift[m] = com[m];
  }
  
  for(i=0; (i<atoms->nr); i++){ 
		 rvec_dec(x0[i],shift);
                 for(m=0;m<3;m++){
		    while(x0[i][m]> box[m][m]/2.) x0[i][m]-=box[m][m];
		    while(x0[i][m]<-box[m][m]/2.) x0[i][m]+=box[m][m];
                 }
  }
}




 

#if 0
void plot_order(real *slDensity[], const char *afile, int nslices,
		  int nr_grps, char *grpname[], real slWidth, 
		  const char **dens_opt,
		  gmx_bool bSymmetrize, const output_env_t oenv,int iterations)
{
  FILE  *den;
  const char *ylabel=NULL;
  int   slice, n,order, incr=(OFF_ORDER2-OFF_ORDER1);
  real  d,dd,density,measure;
  den = fopen("order.xvg", "w");

  //xvgr_legend(den,nr_grps,(const char**)grpname,oenv);
  for (slice = 0; (slice < nslices); slice++) { 
  if(slDensity[RANDOM_PHASE][slice]>0){
    fprintf(den,"%12g  ", slice * slWidth);
    for(order=OFF_ORDER1; order <= OFF_ORDER4; order+=incr){
      for (n = 0; (n < nr_grps); n++) {
        density=slDensity[OFF_NUMBER+n][slice];
        d=( density>0 ? slDensity[order+n][slice]/density : 0.0 ); 
        dd=( density>0 ? (slDensity[order+incr/2+n][slice]/density - d*d)/sqrt(iterations): 0.0 ); 
	fprintf(den,"   %12g %12g", d,dd);
      }
    }
    fprintf(den,"\n");
  }
  }

  gmx_ffclose(den);
}
#endif 



real * load_radii( t_topology *top ) {
	int i;
	real * radii;
        radii = (real * ) malloc( (top->atoms.nr)*sizeof(real )) ;
   	int ntype = top->idef.atnr; 
        for(i=0; (i<top->atoms.nr); i++){
            real vdw=0.,c6,c12,sig6;
            int itype;
            
            itype = top->atoms.atom[i].type;
            c12   = top->idef.iparams[itype*ntype+itype].lj.c12;
            c6    = top->idef.iparams[itype*ntype+itype].lj.c6;
            if ((c6 != 0) && (c12 != 0)) {
            	sig6 = c12/c6;
            	vdw  = pow(sig6,1.0/6.0);
        	}
            radii[i] = 0.5 * vdw;
       }
       return radii;
}

void flip_atoms_in_box(int natoms, int normal,  rvec * x0){
	int i;
	int axis = (normal+1)%3;
	for(i=0;i<natoms;i++){
		x0[i][normal]=-x0[i][normal];
		/* let's preserve chirality */
		x0[i][axis]=-x0[i][axis];
	}
}

#ifdef UNIX
#warning compiling with signal support
void sig_handler(int signo)
{
  if (signo == SIGINT) { 
    printf("\n\n === Caught SIGINT, trying to terminate analysis and save all files...\n\n");
    global_interrupt = 1 ;
  }
}
#endif //UNIX


void calc_intrinsic_density(const char *fn, Flags myflags, atom_id **index, int gnx[],
		  real ***slDensity, int *nslices, int maxlayers, t_topology *top, int ePBC,
		  int axis, int nr_grps, real *slWidth, const output_env_t oenv,
                  real alpha,int *com_opt, const char ** geometry, 
                  int dump_mol,char dens_opt,int ngrps_add){
	enum {NO_ADDITIONAL_INFO=0,ADDITIONAL_INFO=1};
	ITIM * itim;
	real * radii;
  	int curlong=1, totlong=1;
        FILE * surf_cid=NULL;
        FILE * surfmol_cid=NULL;
        FILE * phase_cid=NULL;
        FILE * phase2_cid=NULL;
	int natoms ; 
	real t;
	rvec * x0;
	matrix box;
	int flags;
    	t_trxframe     fr;
        t_trxstatus * status ;	
  	gmx_rmpbc_t  gpbc=NULL;
	t_pbc  pbc;

	radii = load_radii(top);
//        if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
//          gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
#ifdef VIRIAL_EXTENSION
	flags = TRX_NEED_X  |  TRX_READ_V |  TRX_READ_F;
#else
	flags = TRX_NEED_X | TRX_READ_V | TRX_READ_F;
#endif

	if (!read_first_frame(oenv, &status,  fn , &fr, flags)){

	if(dens_opt=='t') { 
          gmx_fatal(FARGS,"Error loading the trajectory: need position, velocities and forces\n");
        } else { 
          gmx_fatal(FARGS,"Error loading the trajectory\n");
	}
}

	for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j] = fr.box[i][j];
        itim = init_intrinsic_surface(myflags,axis, alpha, 0.04,  box, nr_grps, nslices,maxlayers,radii,index,gnx,com_opt,dump_mol,geometry,ngrps_add, top); 

	if(myflags.bDumpPhases){
        	phase_cid=fopen("phase.gro","w");
        	phase2_cid=fopen("phase2.gro","w");
	}


               /* TODO: decide if the density of test lines (0.04) should be hardcoded or not.*/
	FILE * statfile = fopen("stats.dat","w");
	fprintf(statfile, "#column ");
	for(int i = 0 ; i < itim->maxlayers ; i++ ){
		fprintf(statfile, "layer %d surface_density",i);
#ifdef VIRIAL_EXTENSION
		fprintf(statfile, " (layer %d) %d=surface_density %d=dipolez %d=jx jy jz %d=px py pz",i*8,i*8+1,i*8+2,i*8+3,i*8+6);
#endif
	}
	fprintf(statfile, "# \n");

	fclose(statfile);
	if(myflags.bInclusive) check_inclusive_conditions(index, gnx);
  	gpbc = wrap_gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
	do { 
// (SAW) BUG : when no pbc are defined, it loops forever... 
#ifdef UNIX
		if (signal(SIGINT, sig_handler) == SIG_ERR) {
			printf("Error handling SIGINT");
			// sets global_interrupt to 1 when catching SIGINT
		}
#endif //UNIX
		x0 = fr.x;	
        	natoms = fr.natoms;
#if 0
		printf("frame\n");
		for(int i=0;i<natoms;i++) printf("%f %f %f - %f %f %f\n",fr.x[i][0],fr.x[i][1],fr.x[i][2],fr.vir[i][0],fr.vir[i][1],fr.vir[i][2]);
#endif
		for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j] = fr.box[i][j];
    		gmx_rmpbc(gpbc,natoms,box,x0);
        	set_pbc(&pbc, ePBC, box);
                if(myflags.bCenter){
		    put_atoms_in_box(ePBC,box,natoms,x0);
	            /* Make our reference phase whole */
                    if(!myflags.bInclusive) { 
		       remove_phase_pbc(ePBC,&top->atoms,box, x0, axis, index[0], gnx[0]);
		    } else { 
		       remove_phase_pbc(ePBC,&top->atoms,box, x0, axis, index[1], gnx[1]);
	            }
                }
		/* Now center the com of the reference phase in the middle of the box (the one from -box/2 to box/2, 
                   not the gromacs std one:), and shift/rebox the rest accordingly */
                if(!myflags.bInclusive) { 
      			center_coords(&top->atoms,box,x0,axis,index[0],gnx[0]);
		} else { 
      			center_coords(&top->atoms,box,x0,axis,index[1],gnx[1]);
                }
                /* Identify the atom,tops belonging to the intrinsic surface */
	        int success = compute_intrinsic_surface(myflags.bCluster, box, nr_grps, x0, gnx, index,top,natoms,pbc);
		if(itim->info) printf("time=%.3f\n",fr.time);
		if (!success && getenv("NO_IDENTIFICATION_ERROR")==NULL) exit(0);
		if (success) { 
		    if(myflags.bDumpPhases) { 
	                itim->dump_phase_points(INNER_PHASE,top,phase_cid);  // gh
	                itim->dump_phase_points(OUTER_PHASE,top,phase2_cid); // gh
                    }
		    if(myflags.bDump){
	   	        dump_slabs(top,0); 
			//if(bMol)
		    	 //  dump_slabs(top,1);
		    }
		    /* Compute the intrinsic profile */
		       collect_statistics_for_layers(itim,&fr);
	            if(dens_opt!='s') {  
		       //compute_layer_profile(box, index, top, dens_opt, & fr ); 
 		       compute_intrinsic_profile(box, index, top,dens_opt, & fr); 
                    }
		} else {printf("Skipping frame %f\n",fr.time);}

  	} while (read_next_frame(oenv,status,&fr) &&  (global_interrupt == 0) );
        gmx_rmpbc_done(gpbc);

        if(surf_cid  !=NULL)fclose(surf_cid);
        if(surfmol_cid!=NULL)fclose(surfmol_cid);
        if(phase_cid  !=NULL)fclose(phase_cid);
        if(phase2_cid !=NULL)fclose(phase2_cid);

	if(dens_opt!='s')   
	   finalize_intrinsic_profile(slDensity, nslices, slWidth);
#ifdef INTEGER_HISTO
        for(int j=INNER_PHASE;j<itim->nphases;j++){
             char fname[4096];
             sprintf(fname,"intHist.%d.dat",j);
             FILE * cid=fopen(fname,"w");
             dump_int_histo(j,cid);
        }
#endif
	if(dens_opt!='s')free_profile(slDensity);	
}

void free_profile(real *** density){
	int i;
	for(i=0;i<global_itim->n_histo;i++){
	    free ((*density)[i]);
        }
	free(*density);
}

int main(int argc,char *argv[])
{
  const char *desc[] = {
    "Compute partial densities across the box, using an index file. Densities",
    "in kg/m^3, number densities or electron densities can be",
    "calculated. For electron densities, a file describing the number of",
    "electrons for each type of atom should be provided using [TT]-ei[tt].",
    "It should look like:[BR]",
    "   2[BR]",
    "   atomname = nrelectrons[BR]",
    "   atomname = nrelectrons[BR]",
    "The first line contains the number of lines to read from the file.",
    "There should be one line for each unique atom name in your system.",
    "The number of electrons for each atom is modified by its atomic",
    "partial charge."
  };

  output_env_t oenv;
  static real alpha=0.2;
  static const char *dens_opt[] = 
    { NULL, "mass", "number", "charge", "electron", "skip", "tension", "Energy", "U(total energy)",  "x", "y", "z", "X", "Y", "Z","Hcos","Ocos","hcos","ocos","Kinetic",NULL };
  static int  axis = 2;          /* normal to memb. default z  */
  static const char *axtitle="Z"; 
  static const char *geometry[]={NULL,"plane","sphere","cylinder", "generic", NULL}; 
  static int  nslices = 50;      /* nr of slices defined       */
  static int  ngrps   = 1;       /* nr. of groups              */
  static int  ngrps_add  = 0;       /* nr. of groups              */
  static int  maxlayers = 1;       /* max nr. of layers to analyze */
  int com_opt[64];
  int dump_mol=0;

  Flags myflags;
  InitFlags(&myflags);

  t_pargs pa[] = {
    { "-d", FALSE, etSTR, {&axtitle}, 
      "Take the normal on the membrane in direction X, Y or Z." },
    { "-dump", FALSE, etBOOL, {&(myflags.bDump)}, 
      "Dump interfacial atoms." },
    { "-dumpphase", FALSE, etBOOL, {&(myflags.bDumpPhases)}, 
      "Dump phase atoms." },
    { "-sl",  FALSE, etINT, {&nslices},
      "Divide the box in #nr slices." },
    { "-dens",    FALSE, etENUM, {dens_opt},
      "Density"},
    { "-ng",       FALSE, etINT, {&ngrps},
      "Number of groups to compute densities of" },
#ifdef VIRIAL_EXTENSION
    { "-virial",       FALSE, etBOOL, {&(myflags.bVirial)},
      "use virial instead of pressure" },
#endif 
    { "-additional",       FALSE, etINT, {&ngrps_add},
      "Additional groups of which the intrinsic profile per layer should be computed (must belong to the 1st group)" },
    { "-symm",    FALSE, etBOOL, {&(myflags.bSymmetrize)},
      "Symmetrize the density along the axis, with respect to the center. Useful for bilayers." },
    { "-center",  FALSE, etBOOL, {&(myflags.bCenter)},
      "Shift the center of mass along the axis to zero. This means if your axis is Z and your box is bX, bY, bZ, the center of mass will be at bX/2, bY/2, 0."},
    { "-info",  FALSE, etBOOL, {&(myflags.bInfo)},
      "Print additional information during the analysis, such as the number of surface atoms found in each frame."},
    { "-intrinsic",  TRUE, etBOOL, {&(myflags.bIntrinsic)},
      "For back-compatibility only. Does nothing"},
    { "-inclusive",  FALSE, etBOOL, {&(myflags.bInclusive)},
      "Include in the surface calculation molecules from an opposite phase which are not in the main cluster"},
    { "-cluster", FALSE, etBOOL, {&(myflags.bCluster)}, 
      "Filter initial molecules through cluster analysis" },
    { "-alpha", FALSE, etREAL, {&alpha}, 
      "Probe sphere radius for the intrinsic analysis" },
    { "-layers", FALSE, etINT, {&maxlayers}, 
      "number of layers to analyze" },
    { "-MCnorm", FALSE, etBOOL, {&(myflags.bMCnormalization)}, 
      "automatic normalization using MC calculation for arbitrary coordinate systems" },
    { "-mol", FALSE, etBOOL, {&(myflags.bMol)}, 
      "Does molecule-based intrinsic analysis" },
    { "-com", FALSE, etBOOL, {&(myflags.bCom)}, 
      "with the -intrinsic option, perform a molecule-based intrinsic analysis. One should give to this flag the number of atoms in the molecule for each group, space-separated. A zero should be used when no center of mass calculation should be used." },
#if GMX_VERSION >= 50000
    { "-h", FALSE, etBOOL, {&(myflags.bHelp)},  "Print a help message" },
#endif

  };

  real **density;        /* density per slice          */
  real slWidth;          /* width of one slice         */
  char **grpname;        /* groupnames                 */
  char **grpname2;        /* groupnames                 */
  int  nr_electrons;     /* nr. electrons              */
  int  *ngx;             /* sizes of groups            */
  int  *ngx2;             /* sizes of groups            */
  t_electron *el_tab;    /* tabel with nr. of electrons*/
  t_topology *top;       /* topology 		       */ 
  int  ePBC;
  atom_id   **index;     /* indices for all groups     */
  atom_id   **index2;     /* indices for all groups     */
  int  i;


  t_filenm  fnm[] = {    /* files for g_density 	  */
    { efTRX, "-f", NULL,  ffREAD },  
    { efNDX, NULL, NULL,  ffOPTRD }, 
#if GMX_VERSION < 50000
    { efTPX, NULL, NULL,  ffREAD },    	    
#else
    { efTPR, NULL, NULL,  ffREAD },    	    
#endif

    { efDAT, "-ei", "electrons", ffOPTRD }, /* file with nr. of electrons */
    { efXVG,"-o","density",ffWRITE }, 	    
  };

geometry[0]=geometry[1];

#define NFILE asize(fnm)
#if GMX_VERSION >= 50000
  for(int iarg=0; iarg < argc ; iarg ++ ) { // this is horrible, but couldn't get the hang of the new parse_common_args()
	if(!strcmp(argv[iarg],"-h")){
#if defined(NEW_CONTEXT_INTERFACE)
            gmx::CommandLineHelpContext context(&gmx::File::standardError(), gmx::eHelpOutputFormat_Console, NULL,"itim");
#else
            gmx::CommandLineHelpContext context(&gmx::File::standardError(), gmx::eHelpOutputFormat_Console, NULL);
#endif
            gmx::GlobalCommandLineHelpContext global(context);
            gmx_bool res = parse_common_args(&argc,argv,0, NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv) ;

	    exit(0);   
	}
  }
#endif
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME, NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv) ;


  if (myflags.bSymmetrize && !myflags.bCenter) {
    fprintf(stdout,"Can not symmetrize without centering. Turning on -center\n");
    myflags.bCenter = TRUE;
  }
  /* Calculate axis */
  axis = toupper(axtitle[0]) - 'X';
  
#if GMX_VERSION < 50000
  top = read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);     /* read topology file */
#else 
  top = read_top(ftp2fn(efTPR,NFILE,fnm),&ePBC);     /* read topology file */
#endif

  global_masses = (real * ) malloc( (top->atoms.nr)*sizeof(real )) ;
  global_charges = (real * ) malloc( (top->atoms.nr)*sizeof(real )) ;

  for(i=0; (i<top->atoms.nr); i++)  global_masses[i] = top->atoms.atom[i].m;
  for(i=0; (i<top->atoms.nr); i++)  global_charges[i] = top->atoms.atom[i].q;
  if (dens_opt[0][0] != 's'){
     if (dens_opt[0][0] == 'n') {
       for(i=0; (i<top->atoms.nr); i++){
         top->atoms.atom[i].m = 1;  
       }
     } else if (dens_opt[0][0] == 'c') {
       for(i=0; (i<top->atoms.nr); i++)
         top->atoms.atom[i].m = top->atoms.atom[i].q;  
     }
  }

  snew(grpname,ngrps+ngrps_add);
  snew(index,ngrps+ngrps_add);
  snew(ngx,ngrps+ngrps_add);
  snew(CLUSTERCUT,ngrps); 

  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps+ngrps_add,ngx,index,grpname); 
  for(i=0;i<64;i++) com_opt[i]=0;
  if (ngrps>=64) exit(printf("Error, too many groups\n"));
  if(myflags.bCom) { 
     FILE* comfile;
     char comstr[1024],*pchar;
     int iter=0;
     printf("Enter the number of atoms per molecule or residue, or zero if no COM calculation is needed\n");
     for (i=0;i<ngrps;i++){
	printf("group %s: ",grpname[i]);
	scanf("%d",&com_opt[i]);
     }
  } 
  float tmpcut;
  if(myflags.bCluster) {
    for (i=0;i<ngrps;i++) { //gh
      printf("\n Enter a cluster cut-off value for group %s",grpname[i]);
      scanf("%f",&tmpcut);
      CLUSTERCUT[i] = tmpcut;
    } 
  }
  printf("\n");

  if(ngrps==1) { 
	// we assume the system is a one-phase only, or the user wants to caculate the profile 
        // of one phase only: let's play with pointers and fake the two additional groups needed.
	int oldgrps=ngrps;
	ngrps = 3;
        snew(grpname2,ngrps+ngrps_add);
        snew(index2,ngrps+ngrps_add);
        snew(ngx2,ngrps+ngrps_add);
        snew(CLUSTERCUT2,ngrps); 
	for(i=0;i<oldgrps;i++) { 
	   grpname2[i]    = grpname[i];
	   index2[i]      = index[i];
	   ngx2[i]        = ngx[i];
	   CLUSTERCUT2[i] = CLUSTERCUT[i];
	}
	for(i=oldgrps;i<3;i++){
		grpname2[i]=grpname2[i-1];
		index2[i]=index2[i-1];
		CLUSTERCUT2[i]=CLUSTERCUT2[i-1];
		ngx2[i]=ngx2[i-1];
	}
	for(i=3;i<ngrps_add;i++){
		grpname2[i]=grpname[i-3+oldgrps];
		index2[i]=index[i-3+oldgrps];
		CLUSTERCUT2[i]=CLUSTERCUT[i-3+oldgrps];
		ngx2[i]=ngx2[i-3+oldgrps];
	}
	grpname=grpname2;
	index=index2;
	CLUSTERCUT=CLUSTERCUT2;
	ngx=ngx2;
  }

  
  calc_intrinsic_density(ftp2fn(efTRX,NFILE,fnm),myflags,
                           index,ngx,&density,&nslices,maxlayers,
                           top,ePBC,axis,ngrps,&slWidth,oenv,alpha,com_opt,
                           geometry,dump_mol,
                           dens_opt[0][0],ngrps_add);

  plot_intrinsic_density(global_itim->histograms, grpname, opt2fn("-o",NFILE,fnm),dens_opt[0][0]);

  return 0;
}

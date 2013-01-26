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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifdef TIME_PROFILE
#include <sys/time.h>
#endif

#include "sysstuff.h"
#include "string.h"
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "tpxio.h"
#include "physics.h"
#include "gmx_ana.h"
#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

// TODO: make it an option? change the algorithm ? 
#define MAX_KDTREE_CHECK 40


/********************************************************************
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
/* single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk> */



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
#define alloc_resnode()		malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif



struct kdtree *kd_create(int k)
{
	struct kdtree *tree;

	if(!(tree = malloc(sizeof *tree))) {
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
		if(!(node = malloc(sizeof *node))) {
			return -1;
		}
		if(!(node->pos = malloc(dim * sizeof *node->pos))) {
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
			if(!(bptr = buf = malloc(dim * sizeof *bptr))) {
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

	if(!node) return 0;

	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
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
	if(!(rset = malloc(sizeof *rset))) {
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
			if(!(bptr = buf = malloc(dim * sizeof *bptr))) {
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
        
	if(!(rset = malloc(sizeof *rset))) {
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
			if(!(bptr = buf = malloc(dim * sizeof *bptr))) {
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

	if (!(rect = malloc(sizeof(struct kdhyperrect)))) {
		return 0;
	}

	rect->dim = dim;
	if (!(rect->min = malloc(size))) {
		free(rect);
		return 0;
	}
	if (!(rect->max = malloc(size))) {
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
		node = malloc(sizeof *node);
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
#ifdef ALPHA
#include <qhull_tools.h>
#endif
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
#define NORMAL_UNDEFINED -100
typedef	enum { SUPPORT_PHASE=0, INNER_PHASE=1, OUTER_PHASE=2, RANDOM_PHASE=3 } PHASE; // These value are not arbitrary and should not be changed.
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
	real skin;
	real range;
	real *pradii;
	int nphases;
	int ngmxphases;
	int normal;
	int info;
        int n_histo;
	real ** phase;
	int ** phase_index;
	int *n;
	real * alphapoints;
	int nalphapoints;
	int * alpha_index;
#ifdef ALPHA
        setT * alphashape;
#endif
	real * radii;
	int * gmx_alpha_id;
	int * gmx_alpha_phase_id;
	int ** gmx_index;
	int com_opt[64];
        int bCom;
        int bMCnormalization;
        int bOrder;
	real target_mesh_size;
	MESH mesh;
	GEOMETRY geometry;
	PERIODIC *periodic;
	METHOD method;
	char method_name[2][128];
	void (* dump_surface_points)(t_topology*,FILE*);		
	void (* dump_surface_molecules)(t_topology*,FILE*);		
	void (* dump_phase_points)(PHASE,t_topology*,FILE*);		

} ITIM;

ITIM * global_itim;
real * global_real_pointer;
real * global_masses=NULL;

typedef struct {
	real * rdata;
	real size;
	real minsize;
	int nbins;
	int N;
	int iterations;
	real bWidth;
	void (* allocate)(int,int);		
	void (* clear)(int);		
	void (* dump)(int,FILE*);		
	void (* add)(int,real,real);		
	} Histogram;
Histogram * histo_histo;

ITIM * init_intrinsic_surface(int normal, real alpha, real mesh,  matrix box, int ngrps,int *nbins,real * radii, int ** index, int *gnx,int *com_opt, int bOrder, int bMCnormalization, const char ** geometry);

real interpolate_distanceSphere(struct kdtree * Surface, real *P,real *normal, ITIM * itim);
real interpolate_distance3D_2(struct kdtree * Surface, real *P,real *normal, ITIM * itim);
real interpolate_distance3D_3(struct kdtree * Surface, real *P,real *normal, ITIM * itim);

void  compute_intrinsic_profile(matrix box, atom_id **index, t_topology * top);
void  compute_intrinsic_order(matrix box, atom_id **index, t_topology * top);
void  finalize_intrinsic_profile(real *** density, int * nslices, real * slWidth);


void add_histo(int N, real pos, real value){
	int bin ;
	bin = (int) (pos / histo_histo->bWidth);
        if(bin<histo_histo->nbins && bin >=0) histo_histo->rdata[N*histo_histo->nbins+bin] += value;
}
void dump_histo(int N, FILE*file){
	int i;
	for(i=0;i<histo_histo->nbins;i++) fprintf(file,"%f %f\n",i*histo_histo->bWidth,histo_histo->rdata[N*histo_histo->nbins+i]*1.0/histo_histo->iterations);
}
void clear_histo(int N){
	int i;
	for(i=0;i<histo_histo->nbins;i++) histo_histo->rdata[N*histo_histo->nbins+i]=0;
}
void eprintv(real *v, char * s){
	fprintf(stderr,"%f %f %f",v[0],v[1],v[2]);
	fprintf(stderr,"%s",s);
}
void printv(real *v, char * s){
	printf("%f %f %f",v[0],v[1],v[2]);
	printf("%s",s);
}
#ifdef ALPHA
real calculate_circumradius(pointT* p0,pointT* p1,pointT* p2, int dim){
	coordT a = qh_pointdist(p0,p1,dim);
	coordT b = qh_pointdist(p1,p2,dim);
	coordT c = qh_pointdist(p2,p0,dim);

	coordT sum =(a + b + c)*0.5;
	coordT area = sum*(a+b-sum)*(a+c-sum)*(b+c-sum);
	return (real) (a*b*c)/(4*sqrt(area));
}
real calculate_weighted_circumradius(pointT* p,pointT* q,pointT* r, int dim){
	/* this is more complicated than for the sphere, and uses two sqrt...find a better solution */
	
	/* p is the origin, (q-p) points along the new x axis,  (q-p)x(r-p) along the new z azis, 
	   [(q-p)x(r-p)] x (q-p) along the new y axis. */

        realT dx[2],dy[2],dz[2], wq,wp,wr; 
	realT Da,tmp,Dx,Dy;
	real x[3], y[3], z[3],P1[2],P2[2],R1,R2;
	wp = qh_radius(p); wq = qh_radius(q); wr = qh_radius(r);
        wq-=wp; wr-=wp ;
	wp*=wp; wq*=wq; wr*=wr;
        dx[0]=q[0]-p[0]; dx[1]=r[0]-p[0]; 
        dy[0]=q[1]-p[1]; dy[1]=r[1]-p[1]; 
        dz[0]=q[2]-p[2]; dz[1]=r[2]-p[2]; 
	/* let's start by projecting coordinates on the plane of the triangle itself */
	tmp =1./sqrt(dx[0]*dx[0]+dy[0]*dy[0]+dz[0]*dz[0]); 	
	x[0]= dx[0]*tmp;  x[1]= dy[0]*tmp;  x[2]= dz[0]*tmp; 
	z[0]= (x[1]*dz[1] - x[2]*dy[1]); z[1]= (x[2]*dx[1] - x[0]*dz[1]); z[2]= (x[0]*dy[1] - x[1]*dx[1]);
	tmp =1./sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]); 	
	z[0]*=tmp ; z[1]*=tmp ; z[2]*=tmp ;
	/* now x and z are normalized. Let's find y ... */
        y[0]= z[1]*x[2]-z[2]*x[1]; y[1]= z[2]*x[0]-z[0]*x[2]; y[2]= z[0]*x[1]-z[1]*x[0];
	/* now we need only the coordinats in the new (x,y) plane */

	P1[0] =dx[0]*x[0]+dy[0]*x[1]+dz[0]*x[2];
	P1[1] =dx[0]*y[0]+dy[0]*y[1]+dz[0]*y[2];

	P2[0] =dx[1]*x[0]+dy[1]*x[1]+dz[1]*x[2];
	P2[1] =dx[1]*y[0]+dy[1]*y[1]+dz[1]*y[2];

        R1=P1[0]*P1[0]+P1[1]*P1[1]-wq+0*wp;
        R2=P2[0]*P2[0]+P2[1]*P2[1]-wr+0*wp;

	Da = P1[0]*P2[1] - P1[1]*P2[0];
	Dx = - (R1*P2[1] - P1[1]*R2);
	Dy = (R1*P2[0] - P1[0]*R2);
	
	return sqrt((Dx*Dx+Dy*Dy)/(4*(Da*Da)))- sqrt(wp);
}
#endif

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


real interpolate_distanceSphere(struct kdtree * Surface, real *P,real *normal, ITIM * itim){
		              struct kdres * presults;
			      rvec c;
			      real p1[3],p2[3],p3[3],diff[3],dist;
			      int j=0,s2,s1,count=0;
			      /* TODO: comment this: it's the trick about cross prods to determine if a line is within the "cone"*/

			      /* p1 */
			      presults = kd_nearest( Surface, P); 
                              kd_res_item(presults, p1); 
                              rvec_sub(p1,P,diff);
			      kd_res_free(presults);
			      if(norm2(diff) < 1e-6 ) return 0.0;  

			      /* p2 */
			      presults = kd_nearest_range( Surface, p1, 5.*global_itim->alpha); 
			      kd_res_next( presults );
                              kd_res_item(presults, p2); 
                              rvec_sub(p2,P,diff);
                              if(norm2(diff)<1e-6) { kd_res_free(presults); return 0.0 ;  } 

			      cprod(p2,p1,c);
			      s1 = (iprod(c,P)>0?1:-1);
		              while (!kd_res_end(presults)) { 
				   count++;
				   if(count>200){ kd_res_free(presults); return 1e6;}
				   kd_res_next( presults );
                              	   kd_res_item(presults, p3); 
			           cprod(p3,p2,c);
			           s2 = (iprod(c,P)>0?1:-1);
				   if( s2!=s1 ) continue;
			           cprod(p1,p3,c);
			           s2 = (iprod(c,P)>0?1:-1);
				   if( s2!=s1 ) continue;
                                   /* if we reach this point, p3 is the closest particle which has p1 "in the cone" */
                                   rvec_sub(p3,P,diff);
                                   if(norm2(diff)<1e-6) { kd_res_free(presults); return 0.0 ;  } 
				   dist = interpolate_distance3D(p1,p2,p3,P);
		           	   dist = sqrt(iprod(P,P))-dist; 
				   kd_res_free(presults);
				   return dist;
			      }
			      kd_res_free(presults);
			      return 1e6;
}

real interpolate_distance3D(real *A,real *B,real *C,real *I){
/*
let's find r, the vector along I (a given direction) which touches the triangle ABC. 
n = ABxCB ; r=aI ; condition: (r-A).n=0 
-> (aI-A).n=0 -> a=A.n/I.n)  |r|=a|I| = sqrt(I.I) A.n/I.n, if I.n!=0, 
*/
        real tmp;
	int i;
	real ab[3],cb[3],n[3];
        for(i=0;i<3;i++){
		ab[i]=( A[i]-B[i]);
		cb[i]=( C[i]-B[i]);
	}
        cprod(ab,cb,n);
        tmp = iprod(I,n);
	if (fabs(tmp)>1e-8){
//printf("%f %f %f  - ",A[0],A[1],A[2]);
//printf("%f %f %f  - ",B[0],B[1],B[2]);
//printf("%f %f %f  \n ",C[0],C[1],C[2]);
            return (sqrt(iprod(I,I))*iprod(A,n)/tmp); 
        } else { 
	    return sqrt(iprod(A,A));
        }
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

void init_itim_grid(ITIM * itim){
	int i;
	real pos[2];
	for(i=0;i<2;i++){
   	  itim->mesh.n[i] = ceil(itim->box[i]/itim->target_mesh_size);
   	  itim->mesh.size[i] = itim->box[i]/itim->mesh.n[i];
	}
	itim->mesh.nelem = itim->mesh.n[0]*itim->mesh.n[1];
        itim->mesh.flag = (int *) malloc (itim->mesh.nelem * sizeof(int));
	if( itim->mesh.tree != NULL) {
			kd_free( itim->mesh.tree);
			itim->mesh.tree=NULL;
	}
	itim->mesh.tree = kd_create(2);
	pos[0] = - itim->box[0]/2 ; 
	pos[1] = - itim->box[1]/2 ;
	for(i=0;i<itim->mesh.nelem;i++){
             	itim->mesh.flag[i]=DIR_POSITIVE;   
		kd_insert(itim->mesh.tree, pos, (void*)&(itim->mesh.flag[i])); 
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
	itim->alpha_index  = (int *) realloc (itim->alpha_index,(itim->n[0])* sizeof (int)); //SAW TODO: controllare performance qui!
	itim->alphapoints  = (real *) realloc (itim->alphapoints,(itim->n[0])* sizeof (real)*3);
	itim->gmx_alpha_id = (int *) realloc (itim->gmx_alpha_id,(itim->n[0])* sizeof (int));
	partn++;
	presults = kd_nearest_range( itim->mesh.tree, pos, sigma + itim->alpha); 
        while( !kd_res_end( presults ) ) {
          /* get the data and position of the current result item */
          p =  (int*)kd_res_item(presults, res); /* they were init'd all to DIR_POSITIVE. Once we find a 
	  					    touching gridline, we flip the value (see below). */
	  ccc++;
	  if (*p==face) {
		if(flag==0) { 
			flag=1; // if already done, don't add this particle anymore.
		        itim->alpha_index[itim->nalphapoints] = index;
			if(itim->alphapoints==NULL) exit(printf("Error reallocating alphapoints\n"));
//printf("SAW added particle %d to the list of %d. %d lines out of %d done\n",index,itim->nalphapoints+1,counter,itim->mesh.nelem);fflush(stdout);
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
		*p=-1*face;  // flip  the flag of the testline.
		counter ++;
		if(counter == itim->mesh.nelem) {counter = 0 ;  return 0 ;}
	  }
          /* go to the next entry */
          kd_res_next( presults );
        }
//	printf("SAW: face: %d testing particle %d, # of lines in range: %d. So far found %d matches\n",face,partn,ccc,counter);
	kd_res_free( presults );
	return 1;
}


void check_kdtree ( struct kdtree * Surface, real * point , real range, real ** result_points , int * rp_n, int * n_elements ) { // TODO: change the name of this function !!

	struct kdres *presults;
	int iter = 0;
	int offset=*n_elements; // note that this is zero at the beginning, and grows when addressing periodic images.
#ifdef ALPHA
	vertexT* vertex;
#endif
	presults = kd_nearest_range( Surface, point, range); 
 	while( !kd_res_end( presults ) ) {
                  /* get the data and position of the current result item */
		  if(3*(iter+offset)+2 >= *rp_n) { 
			*rp_n+=4*(iter+offset) ; *result_points = realloc (*result_points,*rp_n*sizeof(real) );
			if(*result_points==NULL) exit(printf("Error in realloc\n"));
                  }
#ifdef ALPHA
		  vertex=(vertexT*)kd_res_item( presults, &((*result_points)[3*(iter+offset)]));
                  (*result_points)[3*(iter+offset)+2] = vertex->point[2];
#else
                  (*result_points)[3*(iter+offset)+2] = *((real*)kd_res_item( presults, &((*result_points)[3*(iter+offset)])));
#endif
			// TODO: here goes a realloc, in case ....
		  iter ++ ;
                  /* go to the next entry */
                  kd_res_next( presults );
  		  if(iter> MAX_KDTREE_CHECK) break;
        }
  	kd_res_free( presults );
	*n_elements += iter;
	return ;	
}

int check_kdtree_periodic (struct kdtree * Surface, real * box, real * pos , real range, real ** result_points, int * rp_n ) { 
	real periodic[3];
	int n_elements=0;
	check_kdtree ( Surface, pos, range, result_points ,rp_n, &n_elements) ;
	if(pos[0] < -box[0]/2.  + range  || pos[0] > box[0]/2. - range   ){
		if(pos[0]< -box[0]/2.  +range ){
			periodic[0]=pos[0] + box[0];
		} else { 
			periodic[0]=pos[0] - box[0];
		}
		periodic[1]=pos[1];
		check_kdtree ( Surface, periodic, range, result_points ,rp_n, &n_elements ) ;
	}
	if(pos[1] < -box[1]/2. + range || pos[1] > box[1]/2. - range ){
		if(pos[1]< -box[1]/2. + range){
			periodic[1]=pos[1]+box[1];
		} else { 
			periodic[1]=pos[1]-box[1];
		}
		periodic[0]=pos[0];
		check_kdtree ( Surface, periodic, range, result_points , rp_n , &n_elements) ;
		if(pos[0]< -box[0]/2. + range) 
			periodic[0]=pos[0]+box[0];
		if(pos[0]> box[0]/2. - range )
			periodic[0]=pos[0]-box[0];
		check_kdtree ( Surface, periodic, range, result_points , rp_n,  &n_elements) ;
        }
	return n_elements;
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

	ret *= check_itim_testlines(face,index, pos,sigma,itim,0,gmx_index_phase,top, x0);
	if(pos[0] < border_negative[0] || pos[0] > border_positive[0] ){
		if(pos[0]< border_negative[0]){
			periodic[0]=pos[0]+itim->box[0];
		} else { 
			periodic[0]=pos[0]-itim->box[0];
		}
		periodic[1]=pos[1];
		ret *= check_itim_testlines(face,index,periodic,sigma,itim,0,gmx_index_phase,top,x0);
	}
	if(pos[1] < border_negative[1] || pos[1] > border_positive[1] ) { 
		if(pos[1]< border_negative[1]){
			periodic[1]=pos[1]+itim->box[1];
		} else { 
			periodic[1]=pos[1]-itim->box[1];
		}
		periodic[0]=pos[0];
		ret *= check_itim_testlines(face,index, periodic,sigma,itim,0,gmx_index_phase,top,x0);
		if(pos[0]< border_negative[0])
			periodic[0]=pos[0]+itim->box[0];
		if(pos[0]> border_positive[0])
			periodic[0]=pos[0]-itim->box[0];
		ret *= check_itim_testlines(face,index, periodic,sigma,itim,0,gmx_index_phase,top,x0);
        }
	check_itim_testlines(0,0, NULL,0,NULL,1,NULL,NULL,NULL); // this just resets the flag within check_itim_testlines.
	return ret;
}
int projection(const void * a, const void * b) {
	real a2,b2; //NOTE this is weighted
	int i,j;
	ITIM * itim = global_itim;
	i = *(int*)a;	
	j = *(int*)b;	
	a2 = itim->phase[SUPPORT_PHASE][3*i+2] + itim->radii[i] ;
	b2 = itim->phase[SUPPORT_PHASE][3*j+2] + itim->radii[j] ;
	return (a2>b2?-1:1);
}

// TODO !! NOTE: this allows only computation of surface molecules of the SUPPORT_PHASE
void compute_itim_points(Direction direction, ITIM *itim, int ** gmx_index_phase,t_topology * top, rvec * x0){
	int i=0,index,result;
#ifdef TIME_PROFILE
        struct timeval tp;
        struct timeval tp2;
        gettimeofday(&tp, NULL);
#endif
	qsort((void*)itim->phase_index[SUPPORT_PHASE], itim->n[SUPPORT_PHASE], sizeof(int) , projection);
#ifdef TIME_PROFILE
        gettimeofday(&tp2, NULL);
        fprintf(stderr,"Time to quicksort for itim: millisec=%f\n",1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
        /* The positive side */
        i=0;
	do { 
		index = itim->phase_index[SUPPORT_PHASE][i];
		result = check_itim_testlines_periodic(DIR_POSITIVE,index,&itim->phase[SUPPORT_PHASE][3*index],itim->radii[index],itim,gmx_index_phase,top,x0) ; 
          	i++ ; 
		if(i>=itim->n[SUPPORT_PHASE]) {exit(printf("Error: all (%d) particles scanned, but did not associate all testlines on the positive side...\n",itim->n[SUPPORT_PHASE])); }
	} while (result) ;
	check_itim_testlines(0,0, NULL,0,NULL,1,NULL,NULL,NULL); 
	check_itim_testlines(0,0, NULL,0,NULL,2,NULL,NULL,NULL); 
        /* The negative side */
        i = itim->n[SUPPORT_PHASE]-1;
	do { 
		index = itim->phase_index[SUPPORT_PHASE][i];
		result = check_itim_testlines_periodic(DIR_NEGATIVE,index,&itim->phase[SUPPORT_PHASE][3*index],itim->radii[index],itim,gmx_index_phase,top,x0) ; 
          	i--; 
		if(i<0) {exit(printf("Error: all (%d) particles scanned, but did not associate all testlines on the negative side...\n",itim->n[SUPPORT_PHASE])); }
	} while (result) ;
	check_itim_testlines(0,0, NULL,0,NULL,1,NULL,NULL,NULL); 
	check_itim_testlines(0,0, NULL,0,NULL,2,NULL,NULL,NULL); 

}

#ifdef ALPHA
real determinant3d(realT ** rows) {
 return  determ3d     (rows[0][0], rows[0][1], rows[0][2],
                    rows[1][0], rows[1][1], rows[1][2],
                    rows[2][0], rows[2][1], rows[2][2]);
}


real compute_osculating_sphere_radius(real p[3], real q[3], real r[3] ,real s[3], real wp, real wq, real wr, real ws){
    real M[9],ma[9], det , S[3],D[3],U[3],V[3],tmp,A,B,C,R;      
    real alpha=0, beta=0, gamma=0;
    int i;
    if(ws<0){
	alpha = (q[1]-p[1])*(r[2]-p[2]) - (q[2]-p[2])*(r[1]-p[1]);
	beta  = (q[2]-p[2])*(r[0]-p[0]) - (q[0]-p[0])*(r[2]-p[2]);
	gamma = (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0]);
    }
    for(i=0;i<3;i++){
         ma[i]=p[i]-q[i];
         ma[i+3]=p[i]-r[i];
    }
    if(ws<0){
         ma[6]=alpha; ma[7]=beta; ma[8]=gamma;
    } else {
         for(i=0;i<3;i++)
             ma[i+6]=p[i]-s[i];
    }
    
    tmp = p[0]*p[0]+p[1]*p[1]+p[2]*p[2] - wp*wp;
    S[0]= ( tmp - (q[0]*q[0]+q[1]*q[1]+q[2]*q[2] - wq*wq ))/2;
    S[1]= ( tmp - (r[0]*r[0]+r[1]*r[1]+r[2]*r[2] - wr*wr ))/2;
    if(ws<0) {
	S[2] = alpha * p[0] + beta * p[1] + gamma * p[2] ;
    } else { 
        S[2]= ( tmp - (s[0]*s[0]+s[1]*s[1]+s[2]*s[2] - ws*ws ))/2;
    }
    D[0]= wp-wq; D[1]= wp-wr; if(ws<0){D[2]=0;} else { D[2]= wp-ws;}

    det =   ma[0] * ( ma[4]*ma[8] - ma[7]*ma[5] )
          - ma[1] * ( ma[3]*ma[8] - ma[6]*ma[5] )
          + ma[2] * ( ma[3]*ma[7] - ma[6]*ma[4] );
     M[0] =  ( ma[4]*ma[8] - ma[5]*ma[7] ) / det;
     M[1] = -( ma[1]*ma[8] - ma[7]*ma[2] ) / det;
     M[2] =  ( ma[1]*ma[5] - ma[4]*ma[2] ) / det;
     M[3] = -( ma[3]*ma[8] - ma[5]*ma[6] ) / det;
     M[4] =  ( ma[0]*ma[8] - ma[6]*ma[2] ) / det;
     M[5] = -( ma[0]*ma[5] - ma[3]*ma[2] ) / det;
     M[6] =  ( ma[3]*ma[7] - ma[6]*ma[4] ) / det;
     M[7] = -( ma[0]*ma[7] - ma[6]*ma[1] ) / det;
     M[8] =  ( ma[0]*ma[4] - ma[1]*ma[3] ) / det;
   
     if (fabs(det)<1e-6) return 2.*global_itim->alpha;  /* coplanar / colinear */
 
     for(i=0;i<3;i++){
        V[i] = p[i]- (M[3*i+0]*S[0]+M[3*i+1]*S[1]+M[3*i+2]*S[2]); 
        U[i] = M[3*i+0]*D[0]+M[3*i+1]*D[1]+M[3*i+2]*D[2]; 
     }
     A=(1-U[0]*U[0]-U[1]*U[1]-U[2]*U[2]); 
     B=2*(wp-U[0]*V[0]-U[1]*V[1]-U[2]*V[2]);
     C=wp*wp-V[0]*V[0]-V[1]*V[1]-V[2]*V[2];
     tmp = B*B-4*A*C;
     if (tmp<0) { 
		/* since we eliminated atoms with LJ sigma==0, this means either an
                   atom completely embedded in the other three (i.e. with no exposed surface), 
		   or 4 almost coplanar atoms. Let us return a large (i.e. bigger than alpha) orhtosphere radius in this case */ 
	   return 2.*global_itim->alpha;
     }
     R=(-B+sqrt(B*B-4*A*C))/(2*A);
#if 0
if(!isfinite(R)) exit(printf(":R (%g,%g) not finite: A=%g B=%g C=%g det=%g wp-UV = %g - %g D = %g %g %g disc=%f\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n[%g %g %g %g]\n",
R,tmp,A,B,C,det,wp,U[0]*V[0]-U[1]*V[1]-U[2]*V[2],D[0],D[1],D[2],B*B-4*A*C,p[0],p[1],p[2],q[0],q[1],q[2],r[0],r[1],r[2],s[0],s[1],s[2],wp,wq,wr,ws));
#endif
     tmp=(-B-sqrt(B*B-4*A*C))/(2*A);
     if(R>0){
	if (tmp>0 && tmp < R) R=tmp;
     } else {
	if(tmp<0){
/*
	we assume a particle is within the excluded volume of another one, let's tag it as good for the alpha-complex.
*/
	return tmp;

        }
        R=tmp;
     }
     return R;
}

real compute_osculating_circle_radius(real p[3], real q[3], real r[3] , real wp, real wq, real wr){
	return compute_osculating_sphere_radius(p, q, r , r, wp, wq, wr, -1);
}

real compute_orthosphere_circumradius(real p[3], real q[3], real r[3] ,real s[3], real wp, real wq, real wr, real ws){

  realT *row[3];
  realT dx[3],dy[3],dz[3],dr[3];
  realT Dx,Dy,Dz,Da;
  wq-=wp; wr-=wp ; ws-=wp;
  wp*=wp; wq*=wq; wr*=wr; ws*=ws;
  dx[0]=q[0]-p[0]; dx[1]=r[0]-p[0]; dx[2]=s[0]-p[0];
  dy[0]=q[1]-p[1]; dy[1]=r[1]-p[1]; dy[2]=s[1]-p[1];
  dz[0]=q[2]-p[2]; dz[1]=r[2]-p[2]; dz[2]=s[2]-p[2];
  dr[0]=SQR(dx[0])+SQR(dy[0])+SQR(dz[0])-wq;
  dr[1]=SQR(dx[1])+SQR(dy[1])+SQR(dz[1])-wr;
  dr[2]=SQR(dx[2])+SQR(dy[2])+SQR(dz[2])-ws;

  row[0] = dz; row[1] = dy ; row[2] = dr;  
  Dx = determinant3d(row);
  row[0] = dx; row[1] = dz ; 
  Dy = determinant3d(row);
  row[0] = dy; row[1] = dx ; 
  Dz = determinant3d(row);
  row[0] = dz; row[1] = dy ; row[2] = dx; 
  Da = determinant3d(row);
  if( (Dx*Dx+Dy*Dy+Dz*Dz)/ (4.*Da*Da) < 0 ) exit(printf("Negative value\n")); 
  return ( sqrt( ( (Dx*Dx+Dy*Dy+Dz*Dz)/ (4.*Da*Da) ) ) - sqrt(wp) ) ; 
}

/*
	dim  --> dimension of points
	numpoints --> number of points
	m --> original mesh 
	pm --> new mesh
	alpha --> upper bound for the radius of the empty circumsphere of each simplex 

	compute_alpha_shapes(int dim, int numpoints, MeshModel &m, MeshModel &pm, real alpha, int alphashape)
		build alpha complex or alpha shapes (Edelsbrunner and P.Mucke 1994)from a set of vertices of a mesh with Qhull library (http://www.qhull.org/
		Insert the minimum value of alpha (the circumradius of the triangle) in attribute Quality foreach face.

		The Alpha Shape is the boundary of the alpha complex, that is a subcomplex of the Delaunay triangulation. 
		For a given value of 'alpha', the alpha complex includes all the simplices in the Delaunay 
		triangulation which have an empty circumsphere with radius equal or smaller than 'alpha'. 
		Note that for 'alpha' = 0, the alpha complex consists just of the set P, and for sufficiently large 'alpha', 
		the alpha complex is the Delaunay triangulation DT(P) of P.

		Qhull returns the Delauanay triangulation as a list of tetrahedral facets.

	returns 
		true if no errors occurred;
		false otherwise.
*/
/* NOTE: does not free memory (intentionally) */
setT * compute_alpha_shapes(int dim, int numpoints, coordT* points, real alpha, int *nelements){

        vertexT *vertex,**vertexp;
	

        setT * set, *alphashape=NULL;
        facetT * facet, * neighbor;
        ridgeT *ridge, **ridgep;
        pointT *p0,*p1,*p2;
        real radius;
#ifdef TIME_PROFILE
        struct timeval tp;
        struct timeval tp2;
#endif
	static boolT ismalloc= False;		    /* True if qhull should free points in qh_freeqhull() or reallocation */  
	char flags[]= "qhull d Qt ";	    /* option flags for qhull, see qh_opt.htm   qhull d Qt  */
	int exitcode;			    /* 0 if no error from qhull */
	int convexNumVert,numFacets=0,goodTriangles,vertex_i,vertex_n;
	coordT *center;
		
	int vertexCount=0;
        qh_init_A (stdin, NULL , stderr, 0, NULL);
        qh PROJECTdelaunay=True;
        exitcode= setjmp (qh errexit);
        if (!exitcode){ 
                qh_initflags (flags);
                qh_init_B (points, numpoints, dim, ismalloc); 
#ifdef TIME_PROFILE
        gettimeofday(&tp, NULL);
#endif
                qh_qhull();

#ifdef SUPPORT_IS_SURFACE
	alphashape = qh_settemp(numpoints); 
	FORALLvertices { 
            qh_setappend(&alphashape, vertex); 
        }
        if(alphashape==NULL) exit(printf("Internal error: alphashape not allocated\n"));
	*nelements = numpoints;
        return alphashape;
#endif

#ifdef TIME_PROFILE
        gettimeofday(&tp2, NULL);
        fprintf(stderr,"Time to perform qh_qhull()  millisec=%f\n",1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
                qh_findgood_all (qh facet_list);

	        /* 'qh facet_list' contains the delaunay triangulation */
		//Set facet->center as the Voronoi center
		qh_setvoronoi_all();

		convexNumVert = qh_setsize(qh_facetvertices (qh facet_list, NULL, false));

		//Set of alpha complex triangles for alphashape filtering
		set= qh_settemp(4* (qh num_facets) );  //SAW: this is not deallocated, or?

		qh visit_id++;
		FORALLfacet_(qh facet_list) {
			facet->alphagood=false; 
			facet->seen=true; 
			numFacets++;
			if (!facet->upperdelaunay ) {
				pointT *vertex_points[4]={NULL,NULL,NULL,NULL};
				int vertex_counter=0;
				//For all facets (that are tetrahedrons) calculate the radius of 
                                //the empty circumsphere (orthosphere) considering 
				//the distance between the circumcenter and a vertex of the facet
				FOREACHvertex_(facet->vertices){ 
				     if(vertex_counter==4) { /*we have found a non-simplicial facet, but points should 
                                                               be anyway cocircular: let's take the first 4.*/
						break;
				     }
				     vertex_points[vertex_counter] = vertex->point; 	
				     vertex_counter++;
				}	
				if(vertex_counter<3) exit(printf("Internal error, should not have facets which have less than 3 points"));
                                if(fabs(qh_radius(vertex_points[0])-qh_radius(vertex_points[1]))<1e-6 && 
                                   fabs(qh_radius(vertex_points[1])-qh_radius(vertex_points[2]))<1e-6 &&
                                   fabs(qh_radius(vertex_points[2])-qh_radius(vertex_points[3]))<1e-6) { 
					/* if all 4 atoms have the same radius, then it's easy to compute the 
						   radius of the osculating sphere */
           			        center =  facet->center;
                                        radius =  -qh_radius(vertex_points[0]) + qh_pointdist(vertex_points[0],center,dim);
						
				} else {  /*otherwise we have to solve a linear system ...*/
				        radius =  compute_osculating_sphere_radius(vertex_points[0],vertex_points[1],
					 				           vertex_points[2],vertex_points[3],
				 					           qh_radius(vertex_points[0]),qh_radius(vertex_points[1]),
                                                                                   qh_radius(vertex_points[2]),qh_radius(vertex_points[3]));
				}
	                        if (radius > alpha ) 
				{
					facet->alphagood=false;
					//Compute each ridge (triangle) once and test the cirumference radius with alpha
					facet->visitid= qh visit_id;
#ifndef DO_RIDGES_IN_ALPHA
					qh_makeridges(facet);
					goodTriangles=0;
					FOREACHridge_(facet->ridges) {
						neighbor= otherfacet_(ridge, facet);
						if (( neighbor->visitid != qh visit_id)){ 			
							//Calculate the radius of the circumference 
							p0 = ((vertexT*) (ridge->vertices->e[0].p))->point;
							p1 = ((vertexT*) (ridge->vertices->e[1].p))->point;
							p2 = ((vertexT*) (ridge->vertices->e[2].p))->point;
                                                        if(fabs(qh_radius(p0)-qh_radius(p1))<1e-6 &&  
                                                           fabs(qh_radius(p1)-qh_radius(p2))<1e-6) { 
							        radius = -qh_radius(vertex_points[0])+calculate_circumradius(p0,p1,p2, dim);
                                                        } else  { 
				                                radius =  compute_osculating_circle_radius(p0,p1,p2,
				 					           qh_radius(p0),qh_radius(p1),qh_radius(p2));
                                                        }
							
							if(radius <= alpha){ 
								goodTriangles++;
									//if calculating alpha shape, save the triangle (ridge) for subsequent filtering
							}
						}
					}

					//mark the facet('good' is used as 'marked'). 
					//This facet will have some triangles hidden by the facet's neighbor.
					if( goodTriangles==4) {
						facet->alphagood=true;
							FOREACHridge_(facet->ridges) {
								qh_setappend(&set, ridge); 
				                                FOREACHvertex_i_(ridge->vertices){vertex->seen=vertex->seen2=false;}
							}
					}
#endif
				}
				else //the facet is good. Put all the triangles of the tetrahedron in the mesh
				{
					//Compute each ridge (triangle) once
					facet->visitid= qh visit_id;
					//mark the facet('good' is used as 'marked').
					//This facet will have some triangles hidden by the facet's neighbor.
				        facet->alphagood=true;
					qh_makeridges(facet);
					FOREACHridge_(facet->ridges) {
						neighbor= otherfacet_(ridge, facet);
						if ((neighbor->visitid != qh visit_id)) // let's not process the same ridge twice
                                                        {
								//if calculating alpha shapes, save the triangle for 
                                                                //subsequent filtering	
								qh_setappend(&set, ridge);
								FOREACHvertex_i_(ridge->vertices){
									vertex->seen=false;
									vertex->seen2=false;
                                                                }
							}	
					}
				}
			}
		}
	       // Filter the triangles (only the ones on the boundary of the alpha complex) and build the mesh
			
	        FOREACHridge_(set) { 
			if ( ridge->top->alphagood==false|| 
                             ridge->bottom->alphagood==false ||
                              ridge->top->upperdelaunay || 
                              ridge->bottom->upperdelaunay 
                           ){ 
				FOREACHvertex_i_(ridge->vertices){
                                     if(vertex->seen==false){
					vertexCount++; 
                                        vertex->seen=true;
                                     }
				}
                        }
                }
		alphashape = qh_settemp(vertexCount); 
	        FOREACHridge_(set) {
			if ( ridge->top->alphagood==false|| 
                             ridge->bottom->alphagood==false || 
                              ridge->top->upperdelaunay || 
                              ridge->bottom->upperdelaunay 
                           ){ 
				FOREACHvertex_i_(ridge->vertices){
                                     if(vertex->seen==true && vertex->seen2==false) {
				        qh_setappend(&alphashape, vertex); 
                                        vertex->seen2=true;
                                     }
				}
                        }
                }
	}

//	if (curlong || totlong)
//		fprintf (stderr, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", 
// TODO: fix this!					 totlong, curlong);
        if(alphashape==NULL) exit(printf("Internal error: alphashape not allocated\n"));
	*nelements = vertexCount;
	//if(ismalloc==False) ismalloc=True;
	return alphashape;
}


real interpolate_distance3D_3(struct kdtree * Surface, real *P,real *normal, ITIM * itim){
        real tmp,mindist=1e10,tmpmindist;
	/* NOTE: the alpha parameter here can be used to impose a
	minimum distance between atoms on the surface used to define
	the triangle used to compute the distance. Setting it to
	zero is the same as acceppting all possible atoms */
        real alpha=0;
	int i;
        facetT * facet, ** facetp;
	ridgeT * ridge;
	vertexT* vertex,**vertexp;
	real ba[3],ca[3],da[3],ia[3],n[3],x[3],y[3],Pa[3],pP[2],pa[2],pb[2],pc[2],pi[2],a1,a2,a3,a4,atot;
	real a[3],b[3],diff[3];
	pointT *p[4];
	pointT *pin=NULL;
	struct kdres *pnearest, *presults;
	facetT * neighbor, **neighborp;
	vertexT *Nvertex;
 	real sign=0;
	int near_count=0;
	int isin=0,tmpisin=0,vind=-1,minisin=0;
	/* let's find the nearest point on the surface, and the elements within 2 alpha */
	FORALLfacet_(qh facet_list) {
                        /* let's determine wheter we are inside or outside the surface*/
	    if(facet->alphagood==true){
			realT dist;
			int index=0;
			FOREACHvertex_(facet->vertices){
			    p[index]=vertex->point;
			    index++;
			    if(index>3) break;
                        }
    	        	rvec_sub(p[1],p[0],ba);	
    	        	rvec_sub(p[2],p[0],ca);	
    	        	rvec_sub(p[3],p[0],da);	
    	        	rvec_sub(P,p[0],Pa);	
			cprod(ba,ca,n);
			atot=fabs(iprod(n,da));
			cprod(Pa,ba,n); // Pa x ba .ca 
			a1=fabs(iprod(n,ca));
			cprod(Pa,da,n); // Pa x da .ba 
			a2=fabs(iprod(n,ba));
			cprod(Pa,ca,n); // Pa x ca .da 
			a3=fabs(iprod(n,da));
    	        	rvec_sub(p[2],p[1],ba);	
    	        	rvec_sub(p[3],p[1],ca);	
    	        	rvec_sub(P,p[1],Pa);	
			cprod(Pa,ba,n); // Pa x ba .ca 
			a4=fabs(iprod(n,ca));
			if(fabs(a1+a2+a3+a4-atot)<1e-8) isin=1; // inner side of the surface
            }
        }
	pnearest = kd_nearest_range(Surface,P,1e9);
	i=0;
	while (!kd_res_end(pnearest) ) { 
		int j=0,test[2];
		test[0]=test[1]=1;
       		Nvertex = (vertexT*)kd_res_item(pnearest,NULL); 
                for(j=0;j<i;j++) {
			rvec_sub(Nvertex->point,p[j],diff);	
	        	if (norm2(diff) < alpha*alpha) test[j]=0;
                } 
                if(test[0]*test[1]!=0) { p[i]=Nvertex->point;i++;}
		if(i==3) break;
		kd_res_next( pnearest);
        }
	kd_res_free(pnearest);
	if(i<3) printf("SAW: not found\n");
     	rvec_sub(p[1],p[0],ba);	
    	rvec_sub(p[2],p[0],ca);	
    	cprod(ba,ca,n);
	unitv(n,n);
    	rvec_sub(P,p[0],Pa);	
	tmp = iprod(Pa,n);
	/* let's compute the minimum distance */
        for(i=0;i<3;i++) Pa[i]-=tmp*n[i];
        for(i=0;i<3;i++) x[i]=ba[i];
        unitv(x,x);
        cprod(n,ba,y);
        unitv(y,y);
        pa[0]=0;pa[1]=0;
        pb[0]=iprod(ba,x); pb[1]=iprod(ba,y);
        pc[0]=iprod(ca,x); pc[1]=iprod(ca,y);
        pP[0]=iprod(Pa,x); pP[1]=iprod(Pa,y);
	atot=fabs(triangle_area(pa,pb,pc,global_itim->box));
	a1=  fabs(triangle_area(pa,pb,pP,global_itim->box));
	a2=  fabs(triangle_area(pb,pc,pP,global_itim->box));
	a3=  fabs(triangle_area(pa,pc,pP,global_itim->box));
    	rvec_sub(P,p[0],ba);	
        if(fabs(a1+a2+a3-atot)<1e-8) { 
		tmpmindist=fabs(tmp); 
        } else { 
		tmpmindist=sqrt(iprod(ba,ba));
	}
	if (tmpmindist<mindist){ 
			mindist=tmpmindist;
	}
       

#if oldversion
	FOREACHneighbor_(Nvertex) { 	
	    int index=1;
	    if(neighbor->alphagood==true){
		pin=NULL;
		FOREACHvertex_(neighbor->vertices){
		      if(vertex!=Nvertex){
			if(vertex->seen==true) {
				p[index]=vertex->point; 
				index++;
			} else  {pin=vertex->point;}

                      }
		}
                if(index<3) continue;
    	        rvec_sub(p[1],p[0],ba);	
    	        rvec_sub(p[2],p[0],ca);	
    	        cprod(ba,ca,n);
	        unitv(n,n);
    	        rvec_sub(P,p[0],Pa);	
	        tmp = iprod(Pa,n);
	        /* let's compute the minimum distance */
                for(i=0;i<3;i++) Pa[i]-=tmp*n[i];
                for(i=0;i<3;i++) x[i]=ba[i];
                unitv(x,x);
                cprod(n,ba,y);
                unitv(y,y);
                pa[0]=0;pa[1]=0;
                pb[0]=iprod(ba,x); pb[1]=iprod(ba,y);
                pc[0]=iprod(ca,x); pc[1]=iprod(ca,y);
                pP[0]=iprod(Pa,x); pP[1]=iprod(Pa,y);
	        atot=fabs(triangle_area(pa,pb,pc,global_itim->box));
	        a1=  fabs(triangle_area(pa,pb,pP,global_itim->box));
	        a2=  fabs(triangle_area(pb,pc,pP,global_itim->box));
	        a3=  fabs(triangle_area(pa,pc,pP,global_itim->box));
    	        rvec_sub(P,p[0],ba);	
                if(fabs(a1+a2+a3-atot)<1e-8) { 
	        	tmpmindist=fabs(tmp); 
                } else { 
	        	tmpmindist=sqrt(iprod(ba,ba));
	        }
	        if (tmpmindist<mindist){ 
				mindist=tmpmindist;
		}
            }
	}
#endif // oldversion
#if 0
        if(isin && ddd){
		FILE*cid=fopen("inner.gro","a");
		static int resind=1;
                fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                            resind,"SOL","OW",
			    resind,
			    P[0],P[1],P[2],
                            0.0,0.0,0.0);
		resind++;
		fclose(cid);
	} else {
		FILE*cid=fopen("outer.gro","a");
		static int resind=1;
                fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                            resind,"SOL","OW",
			    resind,
			    P[0],P[1],P[2],
                            0.0,0.0,0.0);
		resind++;
		fclose(cid);
        }
#endif


	if(isin){  return -mindist  ; } else {return mindist;}
}

real interpolate_distance3D_2(struct kdtree * Surface, real *P,real *normal, ITIM * itim){
        real tmp,mindist=1e10,tmpmindist;
	int i;
        facetT * facet, ** facetp;
	ridgeT * ridge;
	vertexT* vertex,**vertexp;
	real ba[3],ca[3],ia[3],n[3],x[3],y[3],Pa[3],pP[2],pa[2],pb[2],pc[2],pi[2],a1,a2,a3,atot;
	real a[3],b[3];
	pointT *p[4];
	pointT *pin=NULL;
	struct kdres *pnearest, *presults;
	facetT * neighbor, **neighborp;
	vertexT *Nvertex;
 	real sign=0;
	int near_count=0;
	int isin=0,tmpisin=0,vind=-1,minisin=0;
	/* let's find the nearest point on the surface, and the elements within 2 alpha */
	presults = kd_nearest_range(Surface,P, 2*global_itim->alpha);
 	while( !kd_res_end( presults ) ) {
		near_count++;
		vertex = (vertexT*)kd_res_item(presults, 0);
		FOREACHneighbor_(vertex) {
			if(neighbor->alphagood==true){
   				realT dist;
    				qh_distplane(P, neighbor, &dist);
				if (dist < qh min_vertex - 2 * qh DISTround) {
        				/* point is clearly inside of facet */
					isin=1;
					break;
    				}

			}
				
 		}
		if(isin){ printf("found in\n"); break;}
		break;
        }
        tmpisin=minisin=isin;
        kd_res_free(presults);
	pnearest = kd_nearest(Surface,P);
	Nvertex = (vertexT*) kd_res_item(pnearest,0);
	tmpisin=isin;
	FOREACHneighbor_(Nvertex) {
	    if(neighbor->alphagood==true){
		vind=0;
		pin=NULL;
                p[0]=Nvertex->point;
                vind++;
		//printf("::: %f %f %f %p %d in:%d\n",Nvertex->point[0],Nvertex->point[1],Nvertex->point[2],(void*)Nvertex->point,Nvertex->seen,isin);
		FOREACHvertex_(neighbor->vertices){
		      if(vertex!=Nvertex){
		    //	printf("%f %f %f %p %d in:%d\n",vertex->point[0],vertex->point[1],vertex->point[2],(void*)vertex->point,vertex->seen,isin);
			if(vertex->seen==true) {
				p[vind]=vertex->point; 
				vind++;
			} else  {pin=vertex->point;}

                      }
		}
                if(vind<3) continue;
    	        rvec_sub(p[1],p[0],ba);	
    	        rvec_sub(p[2],p[0],ca);	
    	        cprod(ba,ca,n);
	        unitv(n,n);
    	        rvec_sub(P,p[0],Pa);	
	        tmp = iprod(Pa,n);
                if(tmpisin==0 && pin !=NULL) {  // otherwise is for sure out
    	        	rvec_sub(pin,p[0],ia);	
    	        	sign=iprod(ia,n);
	        	if((tmp>0)==(sign>0)) tmpisin=1;
#if 0
			if(tmpisin && ddd ){
		    		printf("p0 = %f %f %f \n",0.,0.,0.);
		    		printf("pi = %f %f %f \n",ia[0],ia[1],ia[2]);
		    	 	printf("p1 = %f %f %f \n",ba[0],ba[1],ba[2]);
		    	 	printf("p2 = %f %f %f \n",ca[0],ca[1],ca[2]);
		    		printf("pP = %f %f %f \n",Pa[0],Pa[1],Pa[2]);
		    	 	printf("sign1=%f  sign2=%f\n",sign,tmp);
		    	 	exit(0);
                        }
#endif
                }
	        /* let's compute the minimum distance */
                for(i=0;i<3;i++) Pa[i]-=tmp*n[i];
                for(i=0;i<3;i++) x[i]=ba[i];
                unitv(x,x);
                cprod(n,ba,y);
                unitv(y,y);
                pa[0]=0;pa[1]=0;
                pb[0]=iprod(ba,x); pb[1]=iprod(ba,y);
                pc[0]=iprod(ca,x); pc[1]=iprod(ca,y);
                pP[0]=iprod(Pa,x); pP[1]=iprod(Pa,y);
	        atot=fabs(triangle_area(pa,pb,pc,global_itim->box));
	        a1=  fabs(triangle_area(pa,pb,pP,global_itim->box));
	        a2=  fabs(triangle_area(pb,pc,pP,global_itim->box));
	        a3=  fabs(triangle_area(pa,pc,pP,global_itim->box));
                if(fabs(a1+a2+a3-atot)<1e-8) { 
	        	tmpmindist=sqrt(tmp*tmp); 
                } else { 
	        	tmpmindist=sqrt(iprod(Pa,Pa));
	        }
	        if (tmpmindist<mindist){ 
				mindist=tmpmindist;
				minisin=tmpisin;
		}
            }
	}
        kd_res_free(pnearest);
	if(vind==-1) exit(printf("Error, no facets belonging to the alpha-complex found\n"));
        if(minisin){
		FILE*cid=fopen("inner.gro","a");
		static int resind=1;
                fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                            resind,"SOL","OW",
			    resind,
			    P[0],P[1],P[2],
                            0.0,0.0,0.0);
		resind++;
		fclose(cid);
	} else {
		FILE*cid=fopen("outer.gro","a");
		static int resind=1;
                fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                            resind,"SOL","OW",
			    resind,
			    P[0],P[1],P[2],
                            0.0,0.0,0.0);
		resind++;
		fclose(cid);


        }
	if(tmpisin){  return -mindist  ; } else {return mindist;}
}





#endif


void arrange_alpha_points(ITIM *itim, int ** gmx_index_phase, t_topology * top,rvec * x0){
#ifdef ALPHA
	int nelem;
	setT * alphashape = itim->alphashape;
	vertexT * vertex,**vertexp;
#endif
	int i ;

#ifdef ALPHA
        int nvertices=0;
	
	int nel=0;	
	if(itim->method==METHOD_A_SHAPE) { 
		nelem  = qh_setsize(alphashape);
		if(itim->info)fprintf(stderr,"Number of mesh elements: %d\n",nelem);
		itim->alphapoints = (real *) realloc (itim->alphapoints, nelem * 3 * sizeof(real ));
	        itim->alpha_index  = (int *) realloc (itim->alpha_index, nelem * sizeof (int));
		if(itim->com_opt[SUPPORT_PHASE]==0) { 
		       itim->gmx_alpha_id = (int *) realloc (itim->gmx_alpha_id, nelem * sizeof(int));
		} else {
		       itim->gmx_alpha_id = (int *) realloc (itim->gmx_alpha_id, itim->com_opt[SUPPORT_PHASE] * nelem * sizeof(int));
		}
	}
#endif
	if(Surface!=NULL) {  kd_free(Surface); Surface=NULL; }
	switch(itim->geometry) {
		case SURFACE_PLANE:  Surface  = kd_create( 2 ); break;
#ifdef ALPHA
		case SURFACE_SPHERE : 
		case SURFACE_GENERIC: 
					Surface = kd_create( 3 ); break;
#endif
		default: exit(printf("Geometry not implemented [arrange_alpha_points()]\n"));
	}
	switch(itim->method){
#ifdef ALPHA
	      case METHOD_A_SHAPE:
      	          FOREACHvertex_(alphashape) {
		    nvertices++; 
                        switch(itim->geometry){
		          case SURFACE_PLANE:
		                 if (fabs(vertex->point[0]) > itim->box[0]/2. ||  
                                     fabs(vertex->point[1]) > itim->box[1]/2.) continue ;
                          break;  //SAW CHECK;
		          case SURFACE_SPHERE: 
		          case SURFACE_GENERIC: 
			  break;
		          default: exit(printf("Geometry not implemented [arrange_alpha_points()]\n"));
		              
		        }
      	            	itim->alphapoints[3*nel+0]=vertex->point[0]; 
			itim->alphapoints[3*nel+1]=vertex->point[1]; 
			itim->alphapoints[3*nel+2]=vertex->point[2];
		        if(itim->com_opt[SUPPORT_PHASE] == 0) {
/* TODO check again when using com for the support phase...OW vs Water problem... */
	                          itim->alpha_index[nel]  = qh_id(vertex->point);
			          itim->gmx_alpha_id[nel] = 
                                         gmx_index_phase[SUPPORT_PHASE][ itim->phase_index[SUPPORT_PHASE][itim->alpha_index[nel]]]; 
      	                          kd_insert(Surface, &itim->alphapoints[3*nel], (real*)vertex); //SAW:ref1
//printf("Here: %f %f %f",itim->alphapoints[3*nel],itim->alphapoints[3*nel+1],itim->alphapoints[3*nel+2]);
      	                	  nel++;
			} else {
exit(printf("This has to be rewritten: how to map back to the alpha complex?\n")); /// see SAW:ref1 and SAW:ref2
				int k,ind,dir, pid = itim->phase_index[SUPPORT_PHASE][qh_id(vertex->point)];
				real com[3]={0.0, 0.0, 0.0},mass=0,tot_mass=0,tmpcom[3],firstatom[3];
				/* in case the atoms belong to a molecule we tagged already, skip this point */
                                for(ind=0;ind<nel;ind++){ 
					if(pid == itim->gmx_alpha_id[ind]){
						ind = -1; break;
					}
				}

				for(k=  0 ;  k < itim->com_opt[SUPPORT_PHASE] ; k++) {
                                    /* let's go through all atoms in the tagged molecule ... */
				    ind = itim->com_opt[SUPPORT_PHASE] * ( pid / itim->com_opt[SUPPORT_PHASE] ) + k;
			            /* adding them to our list ...  */		
				    itim->gmx_alpha_id[itim->com_opt[SUPPORT_PHASE]*nel+k] =  gmx_index_phase[SUPPORT_PHASE][ ind ]; 
				    mass = itim->masses[gmx_index_phase[SUPPORT_PHASE][ind]];
				    tot_mass += mass;
                                    /* ...and computing the com. */
				    for(dir=0; dir<3; dir++){ 
					tmpcom[dir] = x0[itim->gmx_alpha_id[ itim->com_opt[SUPPORT_PHASE] * nel + k ]][dir] ;
                                        if(k==0){  
				              firstatom[dir] = tmpcom[dir]; 
                                        } else {
					      while(tmpcom[dir]-firstatom[dir] >  itim->box[dir]/2.) {  tmpcom[dir] -= itim->box[dir]; }
					      while(tmpcom[dir]-firstatom[dir] <  -itim->box[dir]/2.) { tmpcom[dir] += itim->box[dir]; }
	                                }
					com[dir] += tmpcom[dir] * mass;
                                    }
				}
				for(dir=0;dir <3; dir++){
				        com[dir]/=tot_mass;
					while(com[dir] >  itim->box[dir]/2.) com[dir] -= itim->box[dir];
					while(com[dir] < -itim->box[dir]/2.) com[dir] += itim->box[dir];
				        itim->alphapoints[3*nel+dir]=com[dir]; //SAW:ref2
				}
      	                        kd_insert(Surface, &itim->alphapoints[3*nel], &itim->alphapoints[3*nel+2]); 
				nel++;
			 }
/* NOTE:  qh_id() remaps back from the vertex to the index of atoms in the inner phase  (with some periodic copies too!)
	  itim->phase_index remaps the periodic copies to the particles id in the main box.
          gmx_index_phase remaps from the id of the particles  back to the gromacs global atom index.
	  
*/
		  }
                  itim->nalphapoints = nel; 
	          itim->alphapoints = (real *) realloc (itim->alphapoints, itim->nalphapoints * 3 * sizeof(real ));
		  if(itim->com_opt[SUPPORT_PHASE] == 0 ){ 
	  	       itim->gmx_alpha_id = (int *) realloc (itim->gmx_alpha_id, itim->nalphapoints * sizeof(int));
		  } else {  
	  	       itim->gmx_alpha_id = (int *) realloc (itim->gmx_alpha_id, itim->com_opt[SUPPORT_PHASE]*itim->nalphapoints * sizeof(int));
                  }
	      break;
#endif
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
                atom_index=itim->gmx_index[phase][itim->phase_index[phase][i]];
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
void dump_surface_points(t_topology* top,FILE* cid){
	int i,atom_index=0;
        ITIM * itim=global_itim; 
        fprintf(cid,"%s\n",*(top->name));
        fprintf(cid,"%d\n",itim->nalphapoints);
 	for(i=0;i<itim->nalphapoints;i++){
               atom_index = itim->gmx_alpha_id[i];
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
void dump_surface_molecules(t_topology* top,FILE* cid){
	int i,j,atom_index=0,residue_index=0,surface_index=0;
	static int resNR=0;
	static int atomsinRes=0;

        ITIM * itim=global_itim; 
	if(resNR==0) { /* let's count how many residues are in the group. 
			  TODO: maybe this is useful for something else and should be 
			  calculated before and put into the itim structure...*/
		int lastresindex=-1;
		int currentresindex;

		for(i=0;i<itim->n[INNER_PHASE];i++){
			currentresindex=top->atoms.atom[itim->phase_index[INNER_PHASE][i]].resind;
                       // printf("residue: %d %d %s\n",itim->phase_index[INNER_PHASE][i],currentresindex,
	//				*(top->atoms.atomname[itim->phase_index[INNER_PHASE][i]]));
			if(currentresindex!=lastresindex) {
					resNR++;
					lastresindex=currentresindex;
			}
		}
                if(resNR==0){ exit(printf("Error: number of residues in INNER_PHASE is 0\n"));}
		atomsinRes=(itim->n[INNER_PHASE]+1)/resNR;
        }
        if(atomsinRes==0){ exit(printf("Error: number of atoms in INNER_PHASE residues is 0\n"));}
        fprintf(cid,"%s\n",*(top->name));
        fprintf(cid,"%d\n",itim->nalphapoints*atomsinRes); /* SAW This is a pure mess. If two surface atoms are in the same residue, 
							      they will be printed twice. FIXME */
 	for(i=0;i<itim->nalphapoints;i++){
               int offset = (itim->gmx_alpha_id[i]%atomsinRes);
               atom_index = itim->gmx_alpha_id[i] - offset;  /*  here we rewind to the beggining of the residue ... */
               surface_index = itim->phase_index[INNER_PHASE][ itim->alpha_index[i]   ] - offset;
	       for(j=0; j< atomsinRes ; j++) 
                    fprintf(cid,"%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", /* if you gdb up till here, remember that 
									           this broken part of code requires 
									           INNER_PHASE to be the same as SUPPORT_PHASE */
                           top->atoms.atom[atom_index+j].resind+1,
                           *(top->atoms.resinfo[top->atoms.atom[atom_index + j].resind].name),
                           *(top->atoms.atomname[atom_index+j]),
	       	    atom_index+j+1,
	       	    itim->phase[INNER_PHASE][3*(atom_index+j)+0],
	       	    itim->phase[INNER_PHASE][3*(atom_index+j)+1],
	       	    itim->phase[INNER_PHASE][3*(atom_index+j)+2],
                           0.0,0.0,0.0);
        }
        fprintf(cid,"%f %f %f\n",itim->box[0],itim->box[1],itim->box[2]);
}




void init_itim(int nphases) { 
	ITIM * itim;
	int i;
	global_itim=malloc(sizeof(ITIM));
	itim=global_itim;
	sprintf(itim->method_name[METHOD_ITIM],"itim");
	sprintf(itim->method_name[METHOD_A_SHAPE],"alpha-shapes");
        itim->masses=global_masses;
	if(itim->masses==NULL) exit(printf("Internal error, itim->masses not allocated\n"));
	itim->mesh.nelem=0;
	itim->mesh.tree=NULL;
	itim->nphases = nphases;
        itim->ngmxphases=itim->nphases-1;

	itim->n = (int*)malloc(nphases*sizeof(int));
        itim->phase_index =  (int**)malloc(nphases*sizeof(int*));
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
                itim->com_opt[i]=0;
	}
	itim->periodic[INNER_PHASE]=NONE;
	itim->dump_surface_points=dump_surface_points;
	itim->dump_surface_molecules=dump_surface_molecules;
	itim->dump_phase_points=dump_phase_points;
}

void arrange_datapoints( ITIM *itim, rvec * gmx_coords, int *nindex, int ** gmx_index_phase){
/* NOTE:  Here we assume that particles are within a [-box/2:box/2]^3. */
	int atom,additional=0,sign,sign2;
	const real scale = 3.5; /* This defines the thickness of the periodic border, in units of alpha*/
	real radius=0;
	int realloc_factor=2;
	real x,y,z;
	PHASE phase;
        for(phase=SUPPORT_PHASE; phase < itim->ngmxphases ; phase++) { 
          /* at the first iteration, itim->n are all zero  (initialized by init_itim() )*/
          if(itim->n[phase]<nindex[phase]) { 
                itim->n[phase]=nindex[phase];
                itim->phase[phase] = (real * ) realloc(itim->phase[phase],itim->n[phase]*sizeof(real) * 3);
                itim->phase_index[phase] = (int * ) realloc(itim->phase_index[phase],itim->n[phase]*sizeof(int) );
          }
          switch (itim->periodic[phase]) { 
		case NONE:
		case FULL:
          	  switch (itim->geometry) { 
                    case SURFACE_PLANE:
#ifdef ALPHA
                    case SURFACE_SPHERE:
                    case SURFACE_GENERIC:
#endif
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
            		itim->radii=  (real *)realloc(itim->radii, realloc_factor*nindex[phase]*sizeof(real));

                 	for(atom=0 ; atom < nindex[phase] ; atom++) { 
			   x = (gmx_coords)[gmx_index_phase[phase][atom]][(itim->normal+1)%3];
			   y = (gmx_coords)[gmx_index_phase[phase][atom]][(itim->normal+2)%3];
			   z = (gmx_coords)[gmx_index_phase[phase][atom]][itim->normal];
			   radius = itim->radii[atom];

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
                 	   	    itim->radii[nindex[phase]+additional] =  radius; 
			            additional++;
			       }	
			       if(sign*y > itim->box[1]/2.- scale*itim->alpha) {
			            itim->phase[phase][3*(nindex[phase]+additional)] = x;
			            itim->phase[phase][3*(nindex[phase]+additional)+1] = y-sign*itim->box[1];
			            itim->phase[phase][3*(nindex[phase]+additional)+2] = z;
                 	   	    itim->phase_index[phase][nindex[phase]+additional] =  atom; 
                 	   	    itim->radii[nindex[phase]+additional] =  radius; 
			            additional++;
		                    for(sign2=-1; sign2<=1; sign2+=2){ 
			                 if(sign2*x > itim->box[0]/2.- scale*itim->alpha) {
			                      itim->phase[phase][3*(nindex[phase]+additional)] = x-sign2*itim->box[0];
			                      itim->phase[phase][3*(nindex[phase]+additional)+1] = y-sign*itim->box[1];
			                      itim->phase[phase][3*(nindex[phase]+additional)+2] = z;
                 	   	    	      itim->phase_index[phase][nindex[phase]+additional] =  atom; 
                 	   	    	      itim->radii[nindex[phase]+additional] =  radius; 
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
            		                itim->radii=  (real *)realloc(itim->radii, realloc_factor*nindex[phase]*sizeof(real));
			    }

			}
		        itim->n[phase] = nindex[phase] + additional;	
            		itim->radii=  (real *)realloc(itim->radii,itim->n[phase]*sizeof(real));
                        itim->phase[phase] = (real * ) realloc(itim->phase[phase],itim->n[phase]*sizeof(real) * 3);
                        itim->phase_index[phase] = (int * ) realloc(itim->phase_index[phase],itim->n[phase]*sizeof(int) );
			break;
                        default: exit(printf("Geometry not implemented [arrange_datapoints()]\n"));
		  }	
	
		break;
	        default: exit(printf("Periodicity not implemented\n")); break;
	  } // end switch periodicity
#ifdef ALPHA
	  if (itim->method==METHOD_A_SHAPE && phase==SUPPORT_PHASE){ 
          	itim->pradii = qh_set_radii(itim->radii,itim->phase[phase],3,itim->n[phase]);  // SAW !! 
		qh_set_maxval(itim->box);
	  }
#endif
	}
        if(itim->nphases>itim->ngmxphases) itim->n[itim->nphases-1] = itim->n[INNER_PHASE]+itim->n[OUTER_PHASE];
}

Histogram * histo_init(int N, int nbins, real range) { 

	Histogram * histo = histo_histo;
	histo = malloc(N*sizeof(Histogram));
	histo->N=N;
	histo->nbins= nbins; // we are sampling half of the boxlength
	histo->iterations= 0;
	histo->rdata = (real *)calloc(histo->N*histo->nbins,sizeof(real));
        histo->dump = dump_histo;
        histo->clear = clear_histo;
        histo->add = add_histo;
	histo->size = histo->minsize = range;
	histo->bWidth = histo->size / histo->nbins;
	return  histo;
}

static int populate_histogram(int phase, real dist, Histogram * histo, ITIM * itim, real value){
            /* we have values from -box/2 to box/2, let's shift them back to gromacs convention (0:box)*/
	    dist += histo->size/2.;
	    if(dist<histo->size && dist>=0){ 
	 		histo->add(phase,dist,value);
			return 1 ;
	    } 	
	    return 0;
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

real perform_interpolation( struct kdtree *Surface, real * P, real * normal, ITIM * itim) {
		real p1[3],p2[3],p3[3],diff[3];
                rvec v1,v2,vn;
		real a1,a2,a3,atot,dist;
                struct kdres * presults;
		int j=0,found=0;
#ifdef ALPHA
		vertexT* vertex;
#else
		real * zpos;	
#endif

		/* p1 */
		presults = kd_nearest( Surface, P); 
#ifdef ALPHA
                vertex = (vertexT*) kd_res_item(presults, p1); 
		p1[2]= vertex->point[2];
#else 
		zpos    = (real*) kd_res_item(presults, p1); 
		p1[2]=*zpos;
#endif
		kd_res_free(presults);
		presults = kd_nearest_range( Surface, p1, 6.*global_itim->alpha); 
		while(p1[2]*P[2]<0){ /*i.e. they are on opposite sides, this can happen only with ITIM*/
			kd_res_next( presults );
			if(kd_res_end(presults)){  fprintf(stderr,"Warning: no suitable point for interpolation found (if this happens too often something is wrong)\n"); return 1e6 ; }
#ifdef ALPHA
    	                vertex = (vertexT*) kd_res_item(presults, p1); 
			p1[2]= vertex->point[2];
#else 
			zpos    = (real*) kd_res_item(presults, p1); 
			p1[2]=*zpos;
#endif

		}

                rvec_sub(p1,P,diff); if(norm2(diff)<1e-6) { kd_res_free(presults); return 0.0 ;  } 

                do {
			kd_res_next( presults );
			if(kd_res_end(presults)){  fprintf(stderr,"Warning: no suitable point for interpolation found (if this happens too often something is wrong)\n"); return 1e6 ; }
                	kd_res_item(presults, p2); 

#ifdef ALPHA
    	                vertex = (vertexT*) kd_res_item(presults, p2); 
			p2[2]= vertex->point[2];
#else 
			zpos    = (real*) kd_res_item(presults, p2); 
			p2[2]=*zpos;
#endif

		} while(p2[2]*P[2]<0);

                rvec_sub(p2,P,diff); if(norm2(diff)<1e-6) { kd_res_free(presults); return 0.0 ;  } 

		/*Now start iterating over all other (sorted) surface atoms p3, to find a 
		  triangle which p1-p2-p3 encloses our point P, i.e. s.t. a(123)=a(124)+a(234)+a(134)*/
		a1=fabs(triangle_area(p1,p2,P,itim->box));
		while (kd_res_next( presults )) {  
#ifdef ALPHA
    	            vertex = (vertexT*) kd_res_item(presults, p3); 
		    p3[2]= vertex->point[2];
#else 
		    zpos    = (real*) kd_res_item(presults, p3); 
	            p3[2]=*zpos;
#endif


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

void define_support_group(real * radii,int * index, int * gnx, int * com_opt){
	/* excludes atoms with no LJ radius. Used only with alpha shapes*/
        int * tmp,i,j;
	tmp=(int*)malloc((*gnx)*sizeof(int));
	for(i=0,j=0;i<*gnx;i++) {
		if(radii[index[i]]>0.0){
			tmp[j]=index[i];
			j++;
                 }
	}
        *com_opt =  ((*com_opt) * j) / (*gnx)  ; /* a complicated way to say that we remove the 
                                                         number of excluded atom in the molecule. As the 
							 rest of the -com implementation, it assumes that 
							 the index is homogeneous (i.e. same number of atoms 
							 per molecule */
	*gnx=j;
	index=(int*)realloc(index,(*gnx)*sizeof(int));
	for(i=0;i<*gnx;i++) 
		index[i]=tmp[i];
	free(tmp);
}

ITIM * init_intrinsic_surface(int normal, real alpha, real mesh,  matrix box, int ngrps, int *nbins, real * radii, int** gmx_index,int *gnx, int *com_opt,int bOrder, int bMCnormalization , const char ** geometry){
	    ITIM * itim;
	    int i;
            int histofactor=1;
            /* ngrps + 1 here  because of the random phase used to to a MC estimate of the bin volumes for the density profile*/
     	    init_itim(ngrps+1) ;
	    itim = global_itim;
	    itim->info = 1; // TODO put me into cmdline...
#ifndef ALPHA
     	    itim->method = METHOD_ITIM;
#else
	    itim->method = METHOD_A_SHAPE;
#endif
            itim->bOrder=bOrder;
            itim->bMCnormalization=bMCnormalization;
            itim->n_histo=2*itim->nphases;
            /* 1 mass dens x phases  (+random phase),  ( 1 number dens , 4 order, 4 error on order) x phases */
 //TODO: comment this better
            if(itim->bOrder) itim->n_histo = itim->nphases +  (1+4+4) * itim->ngmxphases;
     	    itim->alpha=alpha;
	    itim->skin=0.0;  // TODO: CMDLINE ?
	    itim->range=10.0; // TODO: CMDLINE ? check the other comment about range...
	    switch(geometry[0][0]){
	    	case 'p':  itim->geometry = SURFACE_PLANE; 
     	                   itim->normal=normal; 
      	                   if(itim->method==METHOD_ITIM)  itim->periodic[SUPPORT_PHASE] = NONE; 
#ifdef ALPHA
      	                   if(itim->method==METHOD_A_SHAPE) itim->periodic[SUPPORT_PHASE] = PATCH; 
#endif
			   break;
#ifdef ALPHA
	    	case 's': itim->geometry = SURFACE_SPHERE;
     	                  itim->normal=0; 
                          itim->periodic[SUPPORT_PHASE] = NONE; 
			  break;
	    	case 'g': itim->geometry = SURFACE_GENERIC;
     	                  itim->normal=0; 
                          itim->periodic[SUPPORT_PHASE] = NONE; 
			  break;
#endif
		default:  exit(printf("Geometry not implemented so far.\n"));
	    }

	    itim->target_mesh_size= mesh;
#ifdef ALPHA
            define_support_group(radii,gmx_index[SUPPORT_PHASE],&gnx[SUPPORT_PHASE],&com_opt[SUPPORT_PHASE]);
#endif
            /* itim->com_opt[] are by default zero */
            for(i=0;i<itim->ngmxphases;i++)
	         itim->com_opt[i] = com_opt[i];

	    histo_histo = histo_init(itim->n_histo, *nbins, box[itim->normal][itim->normal] ) ; 
            /* ngrps-1 because we are not computing the density of the SUPPORT_PHASE */
            itim->radii=  (real *)malloc(gnx[SUPPORT_PHASE]*sizeof(real));
	    for(i=0;i<gnx[SUPPORT_PHASE];i++){
		itim->radii[i] = radii[gmx_index[SUPPORT_PHASE][i]];
	    }
            itim->gmx_index = gmx_index;
	    return itim;
}

void finalize_intrinsic_profile(real *** density, int * nslices, real * slWidth){
	int i,j;
	ITIM *itim = global_itim;
	Histogram * histo = histo_histo;
  	*density = (real**) malloc( itim->n_histo* sizeof(real*));
	/* Let's re-determine the number of slices, given the minimum box-size found during the run*/
	*nslices = (int)((histo->minsize*(histo->nbins-1)) / histo->size);  
	*slWidth = histo->bWidth ;
	for(i=0;i<itim->n_histo;i++){
		(*density)[i] = &histo->rdata[ i * histo->nbins];
	}
	for(i=0;i<itim->n_histo;i++){
	   for(j=0;j<histo->nbins;j++){
		   histo->rdata[ i * histo->nbins + j ] = 
                                   (histo->rdata[i* histo->nbins + j]) / 
                                   (histo->iterations * *slWidth);

                   switch(itim->geometry){
			case SURFACE_PLANE:  histo->rdata[ i * histo->nbins + j ]/=(2*(itim->box[0] * itim->box[1]));
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
	     for(i=0;i<itim->n_histo;i++){
                if(i==RANDOM_PHASE) continue;
	   	for(j=0;j<histo->nbins;j++){
                        double norm =  histo->rdata[RANDOM_PHASE*histo->nbins + j ];
                        if(norm>0) histo->rdata[ i * histo->nbins + j ] /= norm ;
                }
             }
        }


printf("box=%f, histo->size= %f = %f, bWidth=%f, nitem=%d =%d, bw*nit=%f\n",itim->box[2],histo->size,histo->minsize,histo->bWidth,histo->nbins,*nslices,histo->nbins*histo->bWidth);
}
void  compute_intrinsic_surface(matrix box, int ngrps, rvec * gmx_coords, int *nindex, atom_id ** gmx_index_phase,t_topology * top){
	
	ITIM * itim = global_itim;
	int ind=itim->normal;
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
	              arrange_datapoints(itim, gmx_coords, nindex,  (int**)gmx_index_phase);
		      compute_itim_points(DIR_POSITIVE,itim,(int**)gmx_index_phase,top,gmx_coords);
		      kd_free(itim->mesh.tree); itim->mesh.tree=NULL;
		break;
#ifdef ALPHA
		case METHOD_A_SHAPE:
	              arrange_datapoints(itim, gmx_coords, nindex,  (int**)gmx_index_phase);
		      itim->alphashape = compute_alpha_shapes(3, itim->n[SUPPORT_PHASE], itim->phase[SUPPORT_PHASE], 
                                                             itim->alpha,&itim->nalphapoints); 
		break;
#endif
		default: exit(printf("Method not implemented yet\n"));
		break;
	}
#ifdef TIME_PROFILE
        gettimeofday(&tp2, NULL);
        fprintf(stderr,"Time to build surface (method: %s): millisec=%f\n",itim->method_name[itim->method],1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
        arrange_alpha_points (itim,(int**)gmx_index_phase,top,gmx_coords); 
	if(itim->info)fprintf(stderr,"Number of surface elements = %d\n",itim->nalphapoints);
}

void compute_layer_profile(matrix box,atom_id ** gmx_index_phase,t_topology * top){  
	int i,j;
	real pos,r;
	ITIM*itim=global_itim;
	Histogram *histo=histo_histo;
        switch (itim->geometry) { 	
	     case SURFACE_PLANE:
        	for(i=0;i<itim->nalphapoints;i++){
        		if(fabs(itim->alphapoints[3*i])<itim->box[0]/2. && fabs(itim->alphapoints[3*i+1])<itim->box[0]/2.){
        			pos=itim->alphapoints[3*i+2];
        			populate_histogram(SUPPORT_PHASE, pos, histo, itim,top->atoms.atom[ itim->gmx_alpha_id[i] ].m);
        	        }
         	
                }
	     break;
             case SURFACE_SPHERE:
             case SURFACE_GENERIC:
                pos=0.0;
         	for(i=0;i<itim->nalphapoints;i++){
			for(j=0;j<3;j++){
                	   r=itim->alphapoints[3*i+j];
			   pos+=r*r;
                        }
			pos=sqrt(pos);
                	populate_histogram(SUPPORT_PHASE, pos, histo, itim,top->atoms.atom[ itim->gmx_alpha_id[i] ].m);
                 }
	     break;
             default: exit(printf("Error: histogram for geometry type %d not implemented\n",itim->geometry)) ;
             break;
        }
}
void compute_histogram(matrix box,atom_id ** gmx_index_phase,t_topology * top){  
/*************************************

          TODO TODO TODO 

This is too cluttered. Reorganize the code...

          TODO TODO TODO 

 ************************************/
    ITIM * itim=global_itim;
    Histogram * histo = histo_histo; 
    static real * result_points = NULL;
    real dist =0;
    int i,j,k,dim;
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
			      size = min(min(itim->box[2],itim->box[1]),itim->box[0])/2; size2=size*size; break;
        default: exit(printf("Not implemented\n"));
    }
// TODO: here gitim takes much more time than itim. Memory issue? check it...
    if(2*size<histo->minsize) histo->minsize=2*size;
    histo->iterations++;
    for(j=INNER_PHASE;j<itim->nphases;j++){
  	/*This is Miguel's original algorithm. */
        if(j==RANDOM_PHASE && !itim->bMCnormalization) continue;
	for(i=0;i<itim->n[j];i++){
                real locmass=0.0,locm,locnumber;
                rvec v1,v2,vm,vn,newpos,oldpos;
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
                            locmass=itim->masses[gmx_index_phase[j][i]];
		            p4[0]=itim->phase[j][3*i]; p4[1]=itim->phase[j][3*i+1]; p4[2]=itim->phase[j][3*i+2];
                     } else { 
                            locmass=itim->box[0]*itim->box[1]*itim->box[2]/itim->n[j];
			    p4[0]=(((real)rand()/RAND_MAX)-0.5)*itim->box[0];
			    p4[1]=(((real)rand()/RAND_MAX)-0.5)*itim->box[1];
			    p4[2]=(((real)rand()/RAND_MAX)-0.5)*itim->box[2];
                     }
                }
                switch(itim->geometry){
                                case SURFACE_PLANE: if(fabs(p4[2])>size) continue; break;
                                case SURFACE_SPHERE: 
                                case SURFACE_GENERIC: 
						    if(iprod(p4,p4)> size2 ) continue; break;
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
		 	       dist = perform_interpolation(Surface, p4, &(normal[0]),itim);
			   break;
#ifdef ALPHA
			   case SURFACE_SPHERE:
			        dist = interpolate_distanceSphere(Surface,p4,&(normal[0]),itim); 
			   break;

			   case SURFACE_GENERIC:
			     	 dist = interpolate_distance3D_3(Surface,p4,&(normal[0]),itim); 
			   break;
#endif // endif ALPHA
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
					case SURFACE_SPHERE :
					case SURFACE_GENERIC:
                                     		v1[0]=p4[0];v1[1]=p4[1];v1[2]=p4[2]; 
           					unitv(v1,v1);
                                     		v2[0]=normal[0];v2[1]=normal[1];v2[2]=normal[2];
           					unitv(v2,v2);
                                              break;
                                        default: exit(printf("Surface %d not implemented\n",itim->geometry)) ; 
                                              break;
                                }
				order[0]=iprod(v1,vm);                                         order2[0]=order[0]*order[0];
				order[1]=iprod(v1,vn); order[1]=(3*order[1]*order[1] -1. )/2.; order2[1]=order[1]*order[1];
				order[2]=iprod(v2,vm);                                         order2[2]=order[2]*order[2];
				order[3]=iprod(v2,vn); order[3]=(3*order[3]*order[3] -1. )/2.; order2[3]=order[3]*order[3];
                 }

                 switch(itim->geometry){ // SAW: TODO: again, this implies -center ... deal with that.
		      //case SURFACE_SPHERE: locnumber=1/(4*M_PI*iprod(p4,p4)); locmass*=locnumber; break;
		      case SURFACE_SPHERE: locnumber=1; break;
		      case SURFACE_GENERIC: locnumber=1; break;
		      case SURFACE_PLANE: locnumber=1; break;
                      default: locnumber=1; break;
                 }
#define NOTEST
#ifndef TEST
		 populate_histogram(j, dist, histo, itim,locmass);
#else 
		 populate_histogram(j, dist, histo, itim,locmass/(4.*M_PI*iprod(p4,p4)));
#endif
               
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
	}
    }
#ifdef TIME_PROFILE
    gettimeofday(&tp2, NULL);
    fprintf(stderr,"Time to make the histogram(s) for %d phases: millisec=%f\n",itim->nphases,1000*(tp2.tv_sec-tp.tv_sec)+((double)tp2.tv_usec-(double)tp.tv_usec)/1000.);
#endif
}

void reset_counters(){
	kd_free(Surface);
	Surface=NULL;
	ITIM * itim = global_itim;
	itim->mesh.nelem=0;
	itim->nalphapoints = 0;
}

void  compute_intrinsic_profile(matrix box, atom_id ** gmx_index_phase, t_topology * top){
	compute_histogram(box,gmx_index_phase,top);
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

  if ( !(in = ffopen(fn,"r")))
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
    (*eltab)[i].atomname = strdup(tempname);
  }
  ffclose(in);
  
  /* sort the list */
  fprintf(stderr,"Sorting list..\n");
  qsort ((void*)*eltab, nr, sizeof(t_electron), 
	 (int(*)(const void*, const void*))compare);

  return nr;
}

void remove_phase_pbc(t_atoms *atoms, matrix box, rvec x0[], int axis, atom_id *index, int index_nr){
     /*we assume here that atoms have been already put into the box */
     /* compute the density at different control points: box edges, middle + some more */
     int rho[5],rho_max,i,nbins=25; 
     static real shift=0.0;
     real bWidth ,z ;
  // TODO: check: this does not work with a solid ... bin must bigger than average interparticle z distance...
     bWidth = box[axis][axis]/nbins;
     while (1){ 
	 if(shift!=0.0){
     	   for(i=0; (i<atoms->nr); i++) {
	   	x0[i][axis]+=shift;
	   }
	   put_atoms_in_box(box,atoms->nr,x0);
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
         if(rho[0] > rho_max/4 || rho[4] > rho_max/4) { /* careful: not >= otherwise the algorithm doesn't work when rho=0 
            					       in all sampled regions */
            	shift+=bWidth*rand()/RAND_MAX;
		fprintf(stderr,"Trying to shift the box by %f nm\n",shift);
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

void calc_electron_density(const char *fn, atom_id **index, int gnx[], 
			   real ***slDensity, int *nslices, t_topology *top,
			   int ePBC,
			   int axis, int nr_grps, real *slWidth, 
			   t_electron eltab[], int nr,gmx_bool bCenter,
                           const output_env_t oenv)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  int natoms;            /* nr. atoms in trj */
  t_trxstatus *status;  
  int i,n,               /* loop indices */
      ax1=0, ax2=0,
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  t_electron *found;     /* found by bsearch */
  t_electron sought;     /* thingie thought by bsearch */
  gmx_rmpbc_t  gpbc=NULL;
 
  real t, 
        z;

  switch(axis) {
  case 0:
    ax1 = 1; ax2 = 2;
    break;
  case 1:
    ax1 = 0; ax2 = 2;
    break;
  case 2:
    ax1 = 0; ax2 = 1;
    break;
  default:
    gmx_fatal(FARGS,"Invalid axes. Terminating\n");
  }

  if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  
  if (! *nslices)
    *nslices = (int)(box[axis][axis] * 10); /* default value */
  fprintf(stderr,"\nDividing the box in %d slices\n",*nslices);

  snew(*slDensity, nr_grps);
  for (i = 0; i < nr_grps; i++)
    snew((*slDensity)[i], *nslices);
  
  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
  /*********** Start processing trajectory ***********/
  do {
    gmx_rmpbc(gpbc,natoms,box,x0);

    if (bCenter)
      center_coords(&top->atoms,box,x0,axis,NULL,0);
    
    *slWidth = box[axis][axis]/(*nslices);
    for (n = 0; n < nr_grps; n++) {      
      for (i = 0; i < gnx[n]; i++) {   /* loop over all atoms in index file */
	  z = x0[index[n][i]][axis];
	  while (z < 0) 
	    z += box[axis][axis];
	  while (z > box[axis][axis])
	    z -= box[axis][axis];
      
	  /* determine which slice atom is in */
	  slice = (z / (*slWidth)); 
	  sought.nr_el = 0;
	  sought.atomname = strdup(*(top->atoms.atomname[index[n][i]]));

	  /* now find the number of electrons. This is not efficient. */
	  found = (t_electron *)
	    bsearch((const void *)&sought,
		    (const void *)eltab, nr, sizeof(t_electron), 
		    (int(*)(const void*, const void*))compare);

	  if (found == NULL)
	    fprintf(stderr,"Couldn't find %s. Add it to the .dat file\n",
		    *(top->atoms.atomname[index[n][i]]));
	  else  
	    (*slDensity)[n][slice] += found->nr_el - 
	                              top->atoms.atom[index[n][i]].q;
	  free(sought.atomname);
	}
    }
      nr_frames++;
  } while (read_next_x(oenv,status,&t,natoms,x0,box));
  gmx_rmpbc_done(gpbc);

  /*********** done with status file **********/
  close_trj(status);
  
/* slDensity now contains the total number of electrons per slice, summed 
   over all frames. Now divide by nr_frames and volume of slice 
*/

  fprintf(stderr,"\nRead %d frames from trajectory. Counting electrons\n",
	  nr_frames);

  for (n =0; n < nr_grps; n++) {
    for (i = 0; i < *nslices; i++)
      (*slDensity)[n][i] = (*slDensity)[n][i] * (*nslices) /
	( nr_frames * box[axis][axis] * box[ax1][ax1] * box[ax2][ax2]);
  }

  sfree(x0);  /* free memory used by coordinate array */
}

void calc_density(const char *fn, atom_id **index, int gnx[], 
		  real ***slDensity, int *nslices, t_topology *top, int ePBC,
		  int axis, int nr_grps, real *slWidth, gmx_bool bCenter,
                  const output_env_t oenv)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  int natoms;            /* nr. atoms in trj */
  t_trxstatus *status;  
  int  **slCount,         /* nr. of atoms in one slice for a group */
      i,j,n,               /* loop indices */
      teller = 0,      
      ax1=0, ax2=0,
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  real t, 
        z;
  char *buf;             /* for tmp. keeping atomname */
  gmx_rmpbc_t  gpbc=NULL;

  switch(axis) {
  case 0:
    ax1 = 1; ax2 = 2;
    break;
  case 1:
    ax1 = 0; ax2 = 2;
    break;
  case 2:
    ax1 = 0; ax2 = 1;
    break;
  default:
    gmx_fatal(FARGS,"Invalid axes. Terminating\n");
  }

  if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  
  if (! *nslices) {
    *nslices = (int)(box[axis][axis] * 10); /* default value */
    fprintf(stderr,"\nDividing the box in %d slices\n",*nslices);
  }
  
  snew(*slDensity, nr_grps);
  for (i = 0; i < nr_grps; i++)
    snew((*slDensity)[i], *nslices);
  
  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
  /*********** Start processing trajectory ***********/
  do {
    gmx_rmpbc(gpbc,natoms,box,x0);

    if (bCenter)
      center_coords(&top->atoms,box,x0,axis,NULL,0);
    
    *slWidth = box[axis][axis]/(*nslices);
    teller++;
    
    for (n = 0; n < nr_grps; n++) {      
      for (i = 0; i < gnx[n]; i++) {   /* loop over all atoms in index file */
	z = x0[index[n][i]][axis];
	while (z < 0) 
	  z += box[axis][axis];
	while (z > box[axis][axis])
	  z -= box[axis][axis];
      
	/* determine which slice atom is in */
	slice = (int)(z / (*slWidth)); 
	(*slDensity)[n][slice] += top->atoms.atom[index[n][i]].m;
      }
    }

    nr_frames++;
  } while (read_next_x(oenv,status,&t,natoms,x0,box));
  gmx_rmpbc_done(gpbc);

  /*********** done with status file **********/
  close_trj(status);
  
  /* slDensity now contains the total mass per slice, summed over all
     frames. Now divide by nr_frames and volume of slice 
     */
  
  fprintf(stderr,"\nRead %d frames from trajectory. Calculating density\n",
	  nr_frames);

  for (n =0; n < nr_grps; n++) {
    for (i = 0; i < *nslices; i++) {
      (*slDensity)[n][i] = (*slDensity)[n][i] * (*nslices) /
	(nr_frames * box[axis][axis] * box[ax1][ax1] * box[ax2][ax2]);
    }
  }

  sfree(x0);  /* free memory used by coordinate array */
}

void plot_density(real *slDensity[], const char *afile, int nslices,
		  int nr_grps, char *grpname[], real slWidth, 
		  const char **dens_opt,
		  gmx_bool bSymmetrize, const output_env_t oenv)
{
  FILE  *den;
  const char *ylabel=NULL;
  int   slice, n, off=0, ncols;
  real  ddd;
  
  ncols=nr_grps;
  switch (dens_opt[0][0]) {
  case 'm': ylabel = "Density (kg m\\S-3\\N)"; break;
  case 'n': ylabel = "Number density (nm\\S-3\\N)"; off=RANDOM_PHASE; ncols=1;  break;
  case 'c': ylabel = "Charge density (e nm\\S-3\\N)"; break;
  case 'e': ylabel = "Electron density (e nm\\S-3\\N)"; break;
  case 's': return;
  }
  
  den = xvgropen(afile, "Partial densities", "Box (nm)", ylabel,oenv);

  xvgr_legend(den,nr_grps,(const char**)grpname,oenv);
  for (slice = 0; (slice < nslices); slice++) { 
    fprintf(den,"%12g  ", slice * slWidth);
    for (n = off; (n < off+ncols); n++) {
      if (bSymmetrize)
	ddd = (slDensity[n][slice]+slDensity[n][nslices-slice-1])*0.5;
      else
	ddd = slDensity[n][slice];
      if (dens_opt[0][0] == 'm')
	fprintf(den,"   %12g", ddd*AMU/(NANO*NANO*NANO));
      else
	fprintf(den,"   %12g", ddd);
    }
    fprintf(den,"\n");
  }

  ffclose(den);
}
 
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

  ffclose(den);
}
 



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
            radii[i] = 0.5 * vdw;  /* SAW !!! */
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

void calc_intrinsic_density(const char *fn, atom_id **index, int gnx[],
		  real ***slDensity, int *nslices, t_topology *top, int ePBC,
		  int axis, int nr_grps, real *slWidth, const output_env_t oenv,
                  real alpha,int *com_opt, int bOrder, const char ** geometry, 
                  int bDump,int bCenter,int bMCnormalization,char dens_opt){
	enum {NO_ADDITIONAL_INFO=0,ADDITIONAL_INFO=1};
	ITIM * itim;
	real * radii;
  	int curlong=1, totlong=1;
        FILE * surf_cid;
        FILE * surfmol_cid;
        FILE * phase_cid;
        FILE * phase2_cid;
	int natoms ; 
	real t;
	rvec * x0;
	matrix box;
        t_trxstatus * status ;	
  	gmx_rmpbc_t  gpbc=NULL;
        surf_cid=fopen("surf.gro","w"); // TODO: from command line, or at least switchable!
        surfmol_cid=fopen("surfmol.gro","w"); 
        phase_cid=fopen("phase.gro","w");
        phase2_cid=fopen("phase2.gro","w");
	radii = load_radii(top);
        if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
          gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  	gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);

        itim = init_intrinsic_surface(axis, alpha, 0.05,  box, nr_grps, nslices,radii,index,gnx,com_opt,bOrder,bMCnormalization,geometry); 
               /* TODO: decide if the density of test lines (0.05) should be hardcoded or not.*/
	do { 
// (SAW) BUG : when no pbc are defined, it loops forever... 
    		gmx_rmpbc(gpbc,natoms,box,x0);
                if(bCenter){
		    put_atoms_in_box(box,natoms,x0);
	            /* Make our reference phase whole */
		    remove_phase_pbc(&top->atoms,box, x0, axis, index[0], gnx[0]);
                }
		/* Now center the com of the reference phase in the middle of the box (the one from -box/2 to box/2, 
                   not the gromacs std one:), and shift/rebox the rest accordingly */
      		center_coords(&top->atoms,box,x0,axis,index[0],gnx[0]);
                /* Identify the atom,tops belonging to the intrinsic surface */
	        compute_intrinsic_surface(box, nr_grps, x0, gnx, index,top);
		if(bDump){
	            itim->dump_phase_points(INNER_PHASE,top,phase_cid); 
	            itim->dump_phase_points(OUTER_PHASE,top,phase2_cid); 
	   	    itim->dump_surface_points(top,surf_cid); 
	   	    itim->dump_surface_molecules(top,surfmol_cid); 
		}
		/* Compute the intrinsic profile */
	        if(dens_opt!='s') {  
		   compute_layer_profile(box, index, top ); 
 		   compute_intrinsic_profile(box, index, top); 
                }
#ifdef ALPHA
		qh_freeqhull(!qh_ALL);  
  		qh_memfreeshort (&curlong, &totlong);
	       	if (curlong || totlong ) if(itim->info)
	    		printf( "internal warning: did not free %d bytes of long memory(%d pieces)\n", totlong, curlong);
#endif


		
  	} while (read_next_x(oenv,status,&t,natoms,x0,box));

	if(dens_opt!='s')   
	   finalize_intrinsic_profile(slDensity, nslices, slWidth);

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
    { NULL, "mass", "number", "charge", "electron", "skip",  NULL };
  static int  axis = 2;          /* normal to memb. default z  */
  static const char *axtitle="Z"; 
  static const char *geometry[]={NULL,"plane","sphere","cylinder", "generic", NULL}; 
  static int  nslices = 50;      /* nr of slices defined       */
  static int  ngrps   = 1;       /* nr. of groups              */
  static gmx_bool bSymmetrize=FALSE;
  int com_opt[64];
  int bCom=0; 
  int bMCnormalization=0; 
  char com_opt_file[1024];
  static gmx_bool bCenter=FALSE;
  static gmx_bool bIntrinsic=FALSE;
  static gmx_bool bDump=FALSE;
  static gmx_bool bOrder=FALSE;
  t_pargs pa[] = {
    { "-d", FALSE, etSTR, {&axtitle}, 
      "Take the normal on the membrane in direction X, Y or Z." },
    { "-dump", FALSE, etBOOL, {&bDump}, 
      "Dump phase(s) and interfacial atoms." },
    { "-sl",  FALSE, etINT, {&nslices},
      "Divide the box in #nr slices." },
    { "-dens",    FALSE, etENUM, {dens_opt},
      "Density"},
    { "-ng",       FALSE, etINT, {&ngrps},
      "Number of groups to compute densities of" },
    { "-symm",    FALSE, etBOOL, {&bSymmetrize},
      "Symmetrize the density along the axis, with respect to the center. Useful for bilayers." },
    { "-center",  FALSE, etBOOL, {&bCenter},
      "Shift the center of mass along the axis to zero. This means if your axis is Z and your box is bX, bY, bZ, the center of mass will be at bX/2, bY/2, 0."},
    { "-intrinsic", FALSE, etBOOL, {&bIntrinsic}, 
      "Perform intrinsic analysis (needs a reference group)" },
    { "-alpha", FALSE, etREAL, {&alpha}, 
      "Probe sphere radius for the intrinsic analysis" },
    { "-MCnorm", FALSE, etBOOL, {&bMCnormalization}, 
      "automatic normalization using MC calculation for arbitrary coordinate systems" },
#ifdef ALPHA
    { "-geometry", FALSE, etENUM, {geometry}, 
      "Geometry of the phase: plane or sphere. Option -d is disregarded for (s)." },
    { "-order", FALSE, etBOOL, {&bOrder}, 
      "Compute order parameter. Switches on -com"},
#endif 
    { "-com", FALSE, etBOOL, {&bCom}, 
      "[only with -alpha option: perform a molecule-based intrinsic analysis. Pass in a string the number of atoms in the molecule of each phase, space-separated. A zero means that no center of mass is used. ]" }
  };

  const char *bugs[] = {
    "When calculating electron densities, atomnames are used instead of types. This is bad.",
  };
  
  real **density;        /* density per slice          */
  real slWidth;          /* width of one slice         */
  char **grpname;        /* groupnames                 */
  int  nr_electrons;     /* nr. electrons              */
  int  *ngx;             /* sizes of groups            */
  t_electron *el_tab;    /* tabel with nr. of electrons*/
  t_topology *top;       /* topology 		       */ 
  int  ePBC;
  atom_id   **index;     /* indices for all groups     */
  int  i;

  t_filenm  fnm[] = {    /* files for g_density 	  */
    { efTRX, "-f", NULL,  ffREAD },  
    { efNDX, NULL, NULL,  ffOPTRD }, 
    { efTPX, NULL, NULL,  ffREAD },    	    
    { efDAT, "-ei", "electrons", ffOPTRD }, /* file with nr. of electrons */
    { efXVG,"-o","density",ffWRITE }, 	    
  };
  
#ifndef ALPHA
geometry[0]=geometry[1];
#endif
#define NFILE asize(fnm)
if(0){
real rad,p[3],q[3],r[3],s[3],wp,wq,wr,ws;
printf("\n-=-=-=-=-=-=-=-=-=-\n");
p[0]=-2.87242; p[1]= -1.21304  ; p[2]=1.66853;
q[0]=-2.85356; q[1]= -0.788321 ; q[2]=1.46017;
r[0]=-2.8576 ; r[1]=-0.571381  ; r[2]=2.01708;
s[0]=-2.88333; s[1]= -1.255    ; s[2]=1.84003;

/*
p[0]= 1.;p[1]= 1.; p[2]= 1.;
q[0]=-1.;q[1]=-1.; q[2]= 1.;
r[0]=-1.;r[1]= 1.; r[2]=-1.;
s[0]= 1.;s[1]=-1.; s[2]=-1.;
*/
//for(wp=0.;wp<sqrt(2);wp+=0.1){
for(wp=0.;wp<0.19;wp+=0.01){
wq=wr=ws=0.;
rad= compute_osculating_sphere_radius(p,q,r,s, wp, wq, wr, ws);
printf("%f %f\n",wp,rad);
}
exit(0);
}

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs,
                    &oenv);

  if (bSymmetrize && !bCenter) {
    fprintf(stderr,"Can not symmetrize without centering. Turning on -center\n");
    bCenter = TRUE;
  }
  /* Calculate axis */
  axis = toupper(axtitle[0]) - 'X';
  
  top = read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);     /* read topology file */


  global_masses = (real * ) malloc( (top->atoms.nr)*sizeof(real )) ;

  for(i=0; (i<top->atoms.nr); i++)  global_masses[i] = top->atoms.atom[i].m;
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
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);
 
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,ngx,index,grpname); 

  if(bCom) { 
     FILE* comfile;
     char comstr[1024],*pchar;
     int iter=0;
     comfile=fopen("masscom.dat","r");
     if(comfile==NULL) exit(printf("Error: file masscom.dat not found\n"));
     fgets(comstr,1024,comfile);
     pchar = strtok(comstr," ") ;
     while (pchar!=NULL ){
           printf("\nhere %s\n",pchar);
           com_opt[iter]=atoi(pchar) ; 
           iter++;
           pchar = strtok(NULL, " ");
     }
  }
  printf("\n%d %d %d\n-----------\n",com_opt[0],com_opt[1],com_opt[2]);
  if (bIntrinsic) { 
    if(ngrps<2) exit(printf("When using -intrinsic please specify at least two groups (can also be the same): the first will be used to compute the intrinsic surface, while the subsequent are used for the density profile calculation.\n"));
    calc_intrinsic_density(ftp2fn(efTRX,NFILE,fnm),index,ngx,&density,&nslices,top,ePBC,axis,ngrps,&slWidth,oenv,alpha,com_opt,bOrder,geometry,bDump,bCenter,bMCnormalization,dens_opt[0][0]);
  } else { 	
    if (dens_opt[0][0] == 'e') {
      nr_electrons =  get_electrons(&el_tab,ftp2fn(efDAT,NFILE,fnm));
      fprintf(stderr,"Read %d atomtypes from datafile\n", nr_electrons);
  
      calc_electron_density(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, 
  			  &nslices, top, ePBC, axis, ngrps, &slWidth, el_tab, 
  			  nr_electrons,bCenter,oenv);
    } else
      calc_density(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, &nslices, top, 
  		 ePBC, axis, ngrps, &slWidth, bCenter,oenv); 
    
  }
    plot_density(density, opt2fn("-o",NFILE,fnm),
	       nslices, ngrps, grpname, slWidth, dens_opt,
	       bSymmetrize,oenv);
#ifdef ALPHA
    if(bOrder)
         plot_order(density, "order.xvg",
	       nslices, ngrps, grpname, slWidth, dens_opt,
	       bSymmetrize,oenv,histo_histo->iterations);
#endif
  do_view(oenv,opt2fn("-o",NFILE,fnm), "-nxy");       /* view xvgr file */
  thanx(stderr);
  return 0;
}

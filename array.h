/*  File: array.h
 *  Authors: Richard Durbin (rd@sanger.ac.uk) and Jean Thierry-Mieg
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1989-
 *  Copyright (C) Genome Research Limited, 1996-
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------
 * This file had its origins as part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description: header for array.c
 * Exported functions:
 *              the Array type and associated macros and functions
 * HISTORY:
 * Last edited: Oct  8 21:53 2014 (rd)
 * * Sep 19 15:42 2014 (rd): switch to long indices to avoid overflow 
 * * May  5 10:57 2013 (rd): new address for rd: rd@sanger.ac.uk
 *-------------------------------------------------------------------
 */

#ifndef ARRAY_DEFINED
#define ARRAY_DEFINED

/* #define ARRAY_CHECK either here or in a single file to
   check the bounds on arr() and arrp() calls
   if defined here can remove from specific C files by defining ARRAY_NO_CHECK
*/

/* #define ARRAY_CHECK */

typedef struct ArrayStruct
  { int   magic ;
    char* base ;    /* char* since need to do pointer arithmetic in bytes */
    long  dim ;     /* length of alloc'ed space */
    int   size ;
    long  max ;     /* largest element accessed via array() */
    int   id ;      /* unique identifier */
  } *Array ;
 
    /* NB we need the full definition for arr() for macros to work
       do not use it in user programs - it is private.
    */

#define ARRAY_MAGIC 8918274

Array   uArrayCreate (long n, int size) ;
#define arrayCreate(n,type)	         uArrayCreate(n,sizeof(type))
Array   uArrayReCreate (Array a, long n, int size) ;
#define arrayReCreate(a,n,type)	         uArrayReCreate(a,n,sizeof(type))
void    arrayDestroy (Array a) ;
Array	arrayCopy (Array a) ;
void    arrayExtend (Array a, long n) ;

     /* array() and arrayp() will extend the array if necessary */

char    *uArray (Array a, long i) ;
#define array(ar,i,type)	(*(type*)uArray(ar,i))
#define arrayp(ar,i,type)	((type*)uArray(ar,i))
char    *uArrayBlock (Array a, long i, long n) ;
#define arrayBlock(ar,i,n,type) ((type*)uArrayBlock(ar,i,n)) /* use when memset(), fread() etc multiple items */

     /* only use arr() when there is no danger of needing expansion */

#if (defined(ARRAY_CHECK) && !defined(ARRAY_NO_CHECK))
char    *uArrCheck (Array a, long i) ;
#define arr(ar,i,type)	(*(type*)uArrCheck(ar,i))
#define arrp(ar,i,type)	((type*)uArrCheck(ar,i))
#else
#define arr(ar,i,type)	((*(type*)((ar)->base + (i)*(ar)->size)))
#define arrp(ar,i,type)	(((type*)((ar)->base + (i)*(ar)->size)))
#endif /* ARRAY_CHECK */

#define arrayMax(ar)  ((ar)->max)

            /* JTM's package to hold sorted arrays of ANY TYPE */
typedef int ArrayOrder(const void*, const void*) ;             /* call back function prototype for sorting arrays */
#define arraySort(a,order)  qsort((a)->base, (a)->max, (a)->size, order)
BOOL    arrayInsert(Array a, void * s, ArrayOrder *order);
BOOL    arrayRemove(Array a, void * s, ArrayOrder *order);
void    arrayCompress(Array a) ;
BOOL    arrayFind(Array a, void *s, long *ip, ArrayOrder *order);

            /* status and memory monitoring */
#define ARRAY_REPORT_MAX 0	/* set to maximum number of arrays to keep track of */
void    arrayStatus (int *nmadep, int* nusedp, 
		     long *memAllocp, long *memUsedp) ; 
/* memUsed only up to REPORT_MAX */
int     arrayReportMark (void) ; /* returns current array number */
void    arrayReport (int j) ;	/* write stderr about all arrays since j */

#endif

/**************************** End of File ******************************/

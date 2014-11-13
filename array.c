/*  File: array.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
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
 * Description:
 *              Arbitrary length arrays
 *              These functions are declared in array.h
 * Exported functions:
 *              See header file: array.h (includes lots of macros)
 * HISTORY:
 * Last edited: Oct  8 21:56 2014 (rd)
 * * Sep 19 15:41 2014 (rd): switch to long indices to avoid overflow 
 * * May  5 10:55 2013 (rd): change RD address to rd@sanger.ac.uk
 * * Feb 14 11:21 2011 (rd): modified in 2009/10 by RD for stand-alone use
 * Created: Thu Dec 12 15:43:25 1989 (mieg)
 *-------------------------------------------------------------------
 */

#include "utils.h"

/********** Array : class to implement variable length arrays **********/

static long totalAllocatedMemory = 0 ;
static int totalNumberCreated = 0 ;
static int totalNumberActive = 0 ;
static Array reportArray = 0 ;

#define arrayExists(a) ((a) && (a)->magic == ARRAY_MAGIC)

Array uArrayCreate (long n, int size)
{
  Array a = mycalloc (1, struct ArrayStruct) ;
  static int isFirst ;

  if (isFirst)
    { isFirst = 0 ;
      if (ARRAY_REPORT_MAX) reportArray = arrayCreate (512, Array) ;
    }
  if (size <= 0) die ("negative size %d in uArrayCreate", size) ;
  if (n < 1)
    n = 1 ;
  totalAllocatedMemory += n * size ;

  a->magic = ARRAY_MAGIC ;
  a->base = _mycalloc (n, size) ;
  a->dim = n ;
  a->max = 0 ;
  a->size = size ;
  a->id = ++totalNumberCreated ;
  ++totalNumberActive ;
  if (reportArray)
    { if (a->id < ARRAY_REPORT_MAX)
	array (reportArray, a->id, Array) = a ;
      else
	{ arrayDestroy (reportArray) ;
	  reportArray = 0 ;
	}
    }
  return a ;
}

/**************/

Array uArrayReCreate (Array a, long n, int size)
{
  if (!arrayExists(a))
    return  uArrayCreate (n, size) ;

  if (a->size != size)
    die ("type size mismatch in arrayReCreate: size %d != a->size %d", size, a->size) ;

  if (n < 1) n = 1 ;

  if (a->dim < n || (a->dim - n)*size > (1 << 20) ) /* free if save > 1 MB */
    { totalAllocatedMemory -= a->dim * size ;
      free (a->base) ;
      a->dim = n ;
      totalAllocatedMemory += a->dim * size ;
      a->base = _mycalloc (n, size) ;
    }
  else
    memset (a->base, 0, n*size) ;

  a->max = 0 ;
  return a ;
}

/**************/

void arrayDestroy (Array a)
{
  if (!arrayExists (a))
    die ("arrayDestroy called on bad array %lx", (long unsigned int) a) ;

  totalAllocatedMemory -= a->dim * a->size ;
  totalNumberActive-- ;
  if (reportArray)
    arr(reportArray, a->id, Array) = 0 ;
  a->magic = 0 ;
  free (a->base) ;
  free (a) ;
}

/**************/

Array arrayCopy (Array a) 
{
  Array new ;

  if (!arrayExists (a)) 
    die ("arrayCopy called on bad array %lx", (long unsigned int) a) ;
 
  new = uArrayCreate (a->dim, a->size) ;
  memcpy (new->base, a->base, a->dim * a->size) ;
  new->max = a->max ;
  return new;
}

/******************************/

void arrayExtend (Array a, long n) 
{
  char *new ;

  if (!arrayExists (a))
    die ("arrayExtend called on bad array %lx", (long unsigned int) a) ;

  if (n < a->dim)
    return ;

  totalAllocatedMemory -= a->dim * a->size ;
  if (a->dim*a->size < 1 << 26)	/* 64MB */
    a->dim *= 2 ;
  else
    a->dim += 1024 + ((1 << 26) / a->size) ;
  if (n >= a->dim)
    a->dim = n + 1 ;

  totalAllocatedMemory += a->dim * a->size ;

  new = _mycalloc (a->dim, a->size) ;
  memcpy (new,a->base,a->size*a->max) ;
  free (a->base) ;
  a->base = new ;

  return;
}

/***************/

char *uArray (Array a, long i)
{
  if (!arrayExists (a))
    die ("array() called on bad array %lx", (long unsigned int) a) ;

  if (i < 0) die("referencing array element %d < 0", i) ;

  if (i >= a->max)
    { if (i >= a->dim)
        arrayExtend (a,i) ;
      a->max = i+1 ;
    }
  return a->base + i*a->size ;
}

/***************/

char *uArrayBlock (Array a, long i, long n)
{
  if (!arrayExists (a))
    die ("arrayBlock() called on bad array %lx", (long unsigned int)a) ;

  if (i < 0) die ("arrayBlock() referencing array element %d < 0", i) ;

  if (i+n >= a->max)
    { if (i+n >= a->dim)
        arrayExtend (a,i+n) ;
      a->max = i+n+1 ;
    }
  return a->base + i*a->size ;
}

/***********************************************/
       /* Finds Entry s from Array  a
        * sorted in ascending order of order()
        * If found, returns TRUE and sets *ip
        * if not, returns FALSE and sets *ip one step left
        */

BOOL arrayFind(Array a, void *s, long *ip, ArrayOrder *order)
{
  int ord ;
  long i = 0 , j, k ;

  if (!arrayExists (a)) 
    die ("arrayFind called on bad array %lx", (long unsigned int) a) ;

  j = arrayMax(a) ;
  if (!j || (ord = order(s,uArray(a,0))) < 0)
    { if (ip) *ip = -1 ;
      return FALSE ;
    }   /* not found */

  if (ord == 0)
    { if (ip) *ip = 0 ;
      return TRUE ;
    }

  if ((ord = order(s,uArray(a,--j))) > 0)
    { if (ip) *ip = j ;
      return FALSE ;
    }
  
  if (ord == 0)
    { if (ip) *ip = j ;
      return TRUE ;
    }

  while(TRUE)
    { k = i + ((j-i) >> 1) ; /* midpoint */
      if ((ord = order(s, uArray(a,k))) == 0)
	{ if (ip) *ip = k ;
	  return TRUE ;
	}
      if (ord > 0) i = k ;
      else j = k ;
      if (i == (j-1))
        break ;
    }
  if (ip)
    *ip = i ;
  return FALSE ;
}

/**************************************************************/
       /* Removes Entry s from Array  a
        * sorted in ascending order of order()
        */

BOOL arrayRemove (Array a, void * s, ArrayOrder *order)
{
  long i;

  if (!arrayExists (a))
    die ("arrayRemove called on bad array %lx", (long unsigned int) a) ;

  if (arrayFind(a, s, &i,order))
    { memmove (uArray(a,i), uArray(a,i+1), (arrayMax(a) - i)*a->size) ;
      arrayMax(a)-- ;
      return TRUE ;
    }
  else
    return FALSE ;
}

/**************************************************************/
       /* Insert Segment s in Array a
        * in ascending order of s.begin
        */

BOOL arrayInsert (Array a, void *s, ArrayOrder *order)
{
  long i ;

  if (!arrayExists (a))
    die ("arrayInsert called on bad array %x", (long unsigned int)a) ;

  if (arrayFind(a, s, &i, order))
    return FALSE ;  /* no doubles */

  uArray (a, arrayMax(a)) ;	/* make space */
  memmove (uArray(a,i+1), uArray(a,i), (arrayMax(a)-1-i)*a->size) ;
  memcpy (uArray(a,i),s,a->size) ;
  return TRUE ;
}

/**************/

void arrayCompress(Array a)
{
  long int i, j, k , as ;
  char *x, *y, *ab ;

  if (!arrayExists (a))
    die ("arrayCompress called on bad array %lx", (long unsigned int) a) ;

  if (arrayMax(a) < 2)
    return ;

  ab = a->base ; 
  as = a->size ;
  for (i = 1, j = 0 ; i < arrayMax(a) ; i++)
    { x = ab + i * as ; y = ab + j * as ;
      for (k = a->size ; k-- ;)		
	if (*x++ != *y++) 
	  goto different ;
      continue ;
      
    different:
      if (i != ++j)
	{ x = ab + i * as ; y = ab + j * as ;
	  for (k = a->size ; k-- ;)	 
	    *y++ = *x++ ;
	}
    }
  arrayMax(a) = j + 1 ;
}

/**************/

int arrayReportMark (void)
{
  if (reportArray)
    return arrayMax (reportArray) ;
  else
    return 0 ;
}

/**************/

void arrayReport (int j)
{
  int i ;
  Array a ;

  fprintf(stderr, "Array report: %d created, %d active, %ld MB allocated\n",   
	  totalNumberCreated, totalNumberActive, totalAllocatedMemory/(1024*1024)) ;

  if (reportArray)
    { i = arrayMax (reportArray) ;
      while (i-- && i > j)
	{ a = arr (reportArray, i, Array) ;
	  if (arrayExists(a))
	    fprintf (stderr, " array %d  size %d, max %ld\n", i, a->size, a->max) ;
	}
    }
}

/**************/

void arrayStatus (int *nmadep, int *nusedp, 
		  long *memAllocp, long *memUsedp)
{ 
  int i ;
  Array a ;

  *nmadep = totalNumberCreated ; 
  *nusedp = totalNumberActive ;
  *memAllocp = totalAllocatedMemory ;
  *memUsedp = 0 ;

  if (reportArray)
    for (i = 0 ; i < arrayMax(reportArray) ; ++i)
      if (arrayExists (a = arr(reportArray, i, Array)))
	*memUsedp += a->max * a->size ;
}

/************************  end of file ********************************/
/**********************************************************************/
 

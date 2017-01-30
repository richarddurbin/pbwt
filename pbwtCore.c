/*  File: pbwtCore.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013-
 *-------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *-------------------------------------------------------------------
 * Description: core functions for pbwt package
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  7 16:21 2015 (rd)
 * * Sep 22 23:02 2014 (rd): change for 64bit arrays
 * Created: Thu Apr  4 11:06:17 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

BOOL isCheck = FALSE ;
BOOL isStats = FALSE ;
DICT *variationDict ;	/* "xxx|yyy" where variation is from xxx to yyy in VCF */

static void pack3init (void) ;	/* forward declaration */

/************** core pbwt options *********************/

void pbwtInit (void) 
{
  variationDict = dictCreate (32) ;
  pack3init () ; 
  sampleInit () ;
  metaDataInit () ;
}

PBWT *pbwtCreate (int M, int N)
{
  PBWT *p = mycalloc (1, PBWT) ; /* cleared so elements default to 0 */

  p->M = M ; p->N = N ;
  p->aFstart = myalloc (M, int) ; int i ; for (i = 0 ; i < M ; ++i) p->aFstart[i] = i ;
  return p ;
}

void pbwtDestroy (PBWT *p)
{
  if (p->chrom) free (p->chrom) ;
  if (p->sites) arrayDestroy (p->sites) ;
  if (p->samples) arrayDestroy (p->samples) ;
  if (p->yz) arrayDestroy (p->yz) ;
  if (p->zz) arrayDestroy (p->zz) ;
  if (p->aFstart) free (p->aFstart) ;
  if (p->aFend) free (p->aFend) ;
  if (p->aRstart) free (p->aRstart) ;
  if (p->aRend) free (p->aRend) ;
  if (p->missingOffset) arrayDestroy (p->missingOffset) ;
  if (p->zMissing) arrayDestroy (p->zMissing) ;
  if (p->dosageOffset) arrayDestroy (p->dosageOffset) ;
  if (p->zDosage) arrayDestroy (p->zDosage) ;
  free (p) ;
}

/*************** subsites, subrange etc. **************/

PBWT *pbwtSubSites (PBWT *pOld, double fmin, double frac)
{
  int M = pOld->M ;
  PBWT *pNew = pbwtCreate (M, 0) ;
  int i, j, k, thresh = M*(1-fmin)  ;
  double bit = 0.0 ;
  uchar *x ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;

  if (!pOld || !pOld->yz) die ("subsites without an existing pbwt") ;
  if (fmin < 0 || fmin >= 1 || frac <= 0 || frac > 1)
    die ("fmin %f, frac %f for subsites out of range\n", fmin, frac) ;

  x = myalloc (M, uchar) ;
  if (pOld->sites) pNew->sites = arrayCreate (4096, Site) ;

  for (i = 0 ; i < pOld->N ; ++i)
    { if ((uOld->c < thresh) && ((bit += frac) > 1.0))
	{ for (j = 0 ; j < M ; ++j) x[uOld->a[j]] = uOld->y[j] ;
	  for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ;
	  pbwtCursorWriteForwards (uNew) ;
	  if (pOld->sites) array(pNew->sites, pNew->N, Site) = arr(pOld->sites, i, Site)  ;
	  ++pNew->N ;
	  bit -= 1.0 ;
	}  
      pbwtCursorForwardsRead (uOld) ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  fprintf (logFile, "subsites with fmin %f, frac %f leaves %d sites\n", fmin, frac, pNew->N) ;

  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pNew->missingOffset = pOld->missingOffset ; pOld->missingOffset = 0 ;
  pNew->zMissing = pOld->zMissing ; pOld->zMissing = 0 ;
  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ;
  return pNew ;
}

PBWT *pbwtSubRange (PBWT *pOld, int start, int end)
{
  int M = pOld->M ;
  PBWT *pNew = pbwtCreate (M, 0) ;
  int i, j, k ;
  uchar *x ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;

  if (!pOld || !pOld->yz) die ("subrange without an existing pbwt") ;
  if (start < 0 || end > pOld->N || end <= start) 
    die ("subrange invalid start %d, end %d", start, end) ;

  x = myalloc (M, uchar) ;
  if (pOld->sites) pNew->sites = arrayCreate (4096, Site) ;

  for (i = 0 ; i < end ; ++i)
    { if (i >= start)
	{ for (j = 0 ; j < M ; ++j) x[uOld->a[j]] = uOld->y[j] ;
	  for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ;
	  pbwtCursorWriteForwards (uNew) ;
	  if (pOld->sites) array(pNew->sites, pNew->N, Site) = arr(pOld->sites, i, Site)  ;
	  ++pNew->N ;
	}  
      pbwtCursorForwardsRead (uOld) ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pNew->missingOffset = pOld->missingOffset ; pOld->missingOffset = 0 ;
  pNew->zMissing = pOld->zMissing ; pOld->zMissing = 0 ;
  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ;
  return pNew ;
}

/***** reverse PBWT - used for local matching to new sequences, and phasing/imputation *****/

void pbwtBuildReverse (PBWT *p)
{
  int i, j, M = p->M ;
  uchar *x = myalloc (M, uchar) ;
  PbwtCursor *uF ;

  if (p->aFend)
    uF = pbwtCursorCreate (p, TRUE, FALSE) ;
  else
    { uF = pbwtCursorCreate (p, TRUE, TRUE) ;
      for (i = 0 ; i < p->N ; ++i)	/* first run forwards to the end */
	pbwtCursorForwardsRead (uF) ;
      pbwtCursorToAFend (uF, p) ;
    }

  /* use p->aFend also to start the reverse cursor - this gives better performance */
  if (!p->aRstart) p->aRstart = myalloc (M, int) ; memcpy (p->aRstart, uF->a, M * sizeof(int)) ;
  p->zz = arrayReCreate (p->zz, arrayMax(p->yz), uchar) ;
  PbwtCursor *uR = pbwtCursorCreate (p, FALSE, TRUE) ; /* will pick up aRstart */
  for (i = p->N ; i-- ; )
    { pbwtCursorReadBackwards (uF) ;
      for (j = 0 ; j < M ; ++j) x[uF->a[j]] = uF->y[j] ;
      for (j = 0 ; j < M ; ++j) uR->y[j] = x[uR->a[j]] ;
      pbwtCursorWriteForwards (uR) ;
    }
  /* save uR->a, which is the lexicographic order of the sequences */
  if (!p->aRend) p->aRend = myalloc (M, int) ; memcpy (p->aRend, uR->a, M * sizeof(int)) ;

  fprintf (logFile, "built reverse PBWT - size %ld\n", arrayMax(p->zz)) ;

  if (isCheck)			/* print out the reversed haplotypes */
    { FILE *fp = fopen ("rev.haps","w") ;
      Array tz = p->yz ; p->yz = p->zz ;
      int* ta = p->aFstart ; p->aFstart = p->aRstart ;
      pbwtWriteHaplotypes (fp, p) ;
      p->yz = tz ; p->aFstart = ta ;
    }

  free (x) ;
  pbwtCursorDestroy (uF) ; pbwtCursorDestroy (uR) ;
}

/*************** make haplotypes ******************/

uchar **pbwtHaplotypes (PBWT *p)	/* NB haplotypes can be space costly */
{
  int M = p->M ;
  int i, j, n = 0 ;
  int *a ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  uchar **hap = myalloc (M, uchar*) ;

  for (i = 0 ; i < M ; ++i) hap[i] = myalloc (p->N, uchar) ;

  for (i = 0 ; i < p->N ; ++i)
    { for (j = 0 ; j < M ; ++j) hap[u->a[j]][i] = u->y[j] ;
      pbwtCursorForwardsRead (u) ;
    }
  pbwtCursorDestroy (u) ;
  return hap ;
}

/****************** Y compression and decompression *************************/
/****************** low level - does not know about PBWT structures *********/

/* pack3 is a three level run length encoding: n times value
   yp & 0x80 = value
   yp & 0x40 == 0 implies n = yp & 0x3f
   yp & 0x40 == 1 implies
     yp & 0x20 == 0 implies n = (yp & 0x1f) << 6
     yp & 0x20 == 1 implies n = (yp & 0x1f) << 11
   This allows coding runs of length up to 64 * 32 * 32 = 64k in 3 bytes.
   Potential factor of ~1000 over a bit array.
   Build a lookup to avoid conditional operation in uncompression.
*/

static int p3decode[128] ;
#define ENCODE_MAX1 64		     /* ~64 */
#define ENCODE_MAX2 ((95-63) << 6)   /* ~1k - is this 32 or 31?*/
#define ENCODE_MAX3 ((127-96) << 11) /* ~64k - ditto */

static void pack3init (void)
{
  int n ;
  for (n = 0 ; n < 64 ; ++n) p3decode[n] = n ;
  for (n = 64 ; n < 96 ; ++n) p3decode[n] = (n-64) << 6 ;
  for (n = 96 ; n < 128 ; ++n) p3decode[n] = (n-96) << 11 ;
}

static inline int pack3Add (uchar yy, uchar *yzp, int n)
/* local utility for pack3 */
{ 
  uchar *yzp0 = yzp ;

  yy <<= 7 ;			/* first move the actual symbol to the top bit */

  while (n >= ENCODE_MAX3) { *yzp++ = yy | 0x7f ; n -= ENCODE_MAX3 ; }
  if (n >= ENCODE_MAX2) { *yzp++ = yy | 0x60 | (n >> 11) ; n &= 0x7ff ; }
  if (n >= ENCODE_MAX1) { *yzp++ = yy | 0x40 | (n >> 6) ; n &= 0x3f ; }
  if (n) { *yzp++ = yy | n ; }
  return yzp - yzp0 ;
}

int pack3 (uchar *yp, int M, uchar *yzp) 
/* pack M chars from y into yz - return number of chars added */
{ 
  int m = 0, m0 ;
  uchar *yzp0 = yzp ;
  uchar ym ;
  
  while (m < M)		/* NB this relies on the M'th character not being 0 or 1 */
    { ym = *yp++ ; m0 = m++ ;
      while (*yp == ym) { ++m ; ++yp ; }
      yzp += pack3Add (ym, yzp, m-m0) ;
    }
  return yzp - yzp0 ;
}

int pack3arrayAdd (uchar *yp, int M, Array ayz)
/* pack M chars onto the end of array ayz - normally use this function */
{
  long max = arrayMax(ayz) ;
  arrayExtend (ayz, max+M) ; /* ensure enough space to copy into */
  int n = pack3 (yp, M, arrp(ayz,max,uchar)) ;
  arrayMax(ayz) = max + n ;
  return n ;
}

int unpack3 (uchar *yzp, int M, uchar *yp, int *n0) 
/* unpack yz into M chars in y - return number of chars unpacked - n0 is number of 0s */
{
  int m = 0 ;
  uchar *yzp0 = yzp ;
  uchar yz ;
  int n ;

  if (n0) *n0 = 0 ;
  while (m < M)
    { yz = *yzp++ ;
      n = p3decode[yz & 0x7f] ;
      m += n ;
      yz = yz >> 7 ;
      if (n0 && !yz) *n0 += n ;
      if (n > 63)		/* only call memset if big enough */
	{ memset (yp, yz, n) ;
	  yp += n ;
	}
      else
	while (n--) *yp++ = yz ;
    }
  if (isCheck && m != M) 
    die ("mismatch m %d != M %d in unpack3 after unpacking %d\n", m, M, yzp-yzp0) ;

  return yzp - yzp0 ;
}

int packCountReverse (uchar *yzp, int M) /* return number of bytes to reverse 1 position */
{ 
  int m = 0 ;
  uchar *yzp0 = yzp ;

  while (m < M)
    m += p3decode[*--yzp & 0x7f] ;
  if (m != M) die ("problem in packCountReverse") ; /* checking assertion */
  return yzp0 - yzp ;
}

#define EATBYTE z = *yzp++ ; n = p3decode[z & 0x7f] ; m += n ; z >>= 7 ; nc[z] += n

int extendMatchForwards (uchar *yzp, int M, uchar x, int *f, int *g)    
		/* x is value to match, [f,g) are start, end of interval */
/* update *f and *g to new match interval, return number of bytes used from yzp */
/* this is more or less standard FM extension, counting nc (=ACC) on fly from BWT */
{
  int m = 0, nc[2] ;
  /* macro for inner loop to move along array */
  uchar z, *yzp0 = yzp ;
  int n = 0 ;

  /* first f */
  nc[0] = nc[1] = 0 ;
  while (m <= *f) { EATBYTE ; }
  if (z == x)
    *f += nc[z] - m ;	        /* equivalent to *f = nc[z] - (m - *f) */
  else
    *f = nc[z] ;

  /* next g */
  if (*g < M)			/* special case g==M so we don't loop beyond M */
    { while (m <= *g) { EATBYTE ; }
      if (z == x)
	*g += nc[z] - m ;
      else
	*g = nc[z] ;
    }

  while (m < M) { EATBYTE ; }	/* complete reading the column */
  
  if (*g == M)          	/* here is the special case - must come after end of column */
    *g = x ? (M-nc[0]) : nc[0] ;

  if (x)			/* add on nc[0] because in block of 1s which follows 0s */
    { *f += nc[0] ;
      *g += nc[0] ;
    }

  return yzp - yzp0 ;
}

int extendPackedForwards (uchar *yzp, int M, int *f, uchar *zp)
{
  int m = 0, nc[2], n ;
  uchar z, *yzp0 = yzp ;

  nc[0] = nc[1] = 0 ;
  while (m <= *f) { EATBYTE ; }	/* find the block containing *f */
  *f += nc[z] - m ;	        /* equivalent to *f = nc[z] - (m - *f) */
  *zp = z ;			/* save the value */
  while (m < M) { EATBYTE ; }	/* complete the column */
  if (*zp) *f += nc[0] ;	/* add on c because in block of 1s which follows 0s */

  return yzp - yzp0 ;
}

int extendPackedBackwards (uchar *yzp, int M, int *f, int c, uchar *zp)
{
  int n, m = 0, nc[2] ;
  uchar z, *yzp0 = yzp, *yzp1 ; 

  while (m < M)			/* first go back to start of previous block */
    m += p3decode[*--yzp & 0x7f] ;
  yzp1 = yzp ;			/* record the start so we can return the difference */

  m = 0 ; nc[0] = nc[1] = 0 ;
  if (*f < c)			/* it was a 0 */
    { while (nc[0] <= *f) { EATBYTE ; }
      *f += nc[1] ;	        /* equivalent to *f = m - (nc[0] - *f) */
      *zp = 0 ;
    }
  else				/* it was a 1 */
    { while (nc[1] <= *f-c) { EATBYTE ; }
      *f += nc[0] - c ;	        /* equivalent to *f = m - (nc[1] - (*f-c)) */
      *zp = 1 ;
    }
  				/* we don't need to finish the block this time */
  return yzp0 - yzp1 ;
}

/************ block extension algorithms, updating whole arrays ************/
/* we could do these also on the packed array with memcpy */

PbwtCursor *pbwtNakedCursorCreate (int M, int *aInit) 
{
  PbwtCursor *u = mycalloc (1, PbwtCursor) ;
  int i ; 

  u->M = M ;
  u->y = myalloc (M+1, uchar) ; u->y[M] = Y_SENTINEL ;
  u->a = myalloc (M, int) ;
  if (aInit) memcpy (u->a, aInit, M*sizeof(int)) ;
  else for (i = 0 ; i < M ; ++i) u->a[i] = i ;
  u->b = myalloc (M, int) ;
  u->d = mycalloc (M+1, int) ;
  u->d[0] = 1 ; u->d[u->M] = 1 ; /* sentinel values */
  u->e = myalloc (M+1, int) ;
  u->u = myalloc (M+1, int) ;
  return u ;
}

PbwtCursor *pbwtCursorCreate (PBWT *p, BOOL isForwards, BOOL isStart)
{
  if (isForwards && !p->yz) p->yz = arrayCreate (1<<20, uchar) ;
  if (!isForwards && !p->zz) p->zz = arrayCreate (1<<20, uchar) ;
  PbwtCursor *u ;
  if (isForwards && isStart) u = pbwtNakedCursorCreate (p->M, p->aFstart) ; 
  else if (isForwards && !isStart) u = pbwtNakedCursorCreate (p->M, p->aFend) ; 
  else if (!isForwards && isStart) u = pbwtNakedCursorCreate (p->M, p->aRstart) ;
  else if (!isForwards && !isStart) u = pbwtNakedCursorCreate (p->M, p->aRend) ;
  if (isForwards) u->z = p->yz ; else u->z = p->zz ;
  if (isStart) 
    if (arrayMax(u->z))
      { u->nBlockStart = 0 ;
	u->n = unpack3 (arrp(u->z,0,uchar), p->M, u->y, &u->c) ;
	u->isBlockEnd = TRUE ;
      }
    else 
      { u->n = 0 ; 
	u->isBlockEnd = FALSE ;
      }
  else 				/* isEnd */
    { u->n = arrayMax(u->z) ;
      u->isBlockEnd = FALSE ;
    }
  return u ;
}

void pbwtCursorDestroy (PbwtCursor *u)
{
  free (u->y) ;
  free (u->a) ;
  free (u->b) ;
  free (u->d) ;
  free (u->e) ;
  free (u->u) ;
  free (u) ;
}

void pbwtCursorForwardsA (PbwtCursor *x) /* algorithm 1 in the manuscript */
{
  int u = 0, v = 0 ;
  int i ;
    
  for (i = 0 ; i < x->M ; ++i)
    if (x->y[i] == 0)
      x->a[u++] = x->a[i] ;
    else			/* y[i] == 1, since bi-allelic */
      x->b[v++] = x->a[i] ;

  memcpy (x->a+u, x->b, v*sizeof(int)) ;
}

void pbwtCursorBackwardsA (PbwtCursor *x) /* undo algorithm 1 */
{
  int u = 0, v = 0 ;
  int i ;
  int *t = x->a ; x->a = x->b ; x->b = t ; /* will copy back from b to a */

  for (i = 0 ; i < x->M ; ++i)
    if (x->y[i] == 0)
      x->a[i] = x->b[u++] ;
    else			/* y[i] == 1, since bi-allelic */
      x->a[i] = x->b[x->c + v++] ;
}

void pbwtCursorForwardsAD (PbwtCursor *x, int k) /* algorithm 2 in the manuscript */
{
  int u = 0, v = 0 ;
  int i ;
  int p = k+1 ; int q = k+1 ;
    
  for (i = 0 ; i < x->M ; ++i)
    { if (x->d[i] > p) p = x->d[i] ;
      if (x->d[i] > q) q = x->d[i] ;
      if (x->y[i] == 0)		/* NB x[a[i]] = y[i] in manuscript */
	{ x->a[u] = x->a[i] ;
	  x->d[u] = p ;
	  ++u ; p = 0 ;
	}
      else			/* y[i] == 1, since bi-allelic */
	{ x->b[v] = x->a[i] ;
	  x->e[v] = q ;
	  ++v ; q = 0 ;
	}
    }

  memcpy (x->a+u, x->b, v*sizeof(int)) ;
  memcpy (x->d+u, x->e, v*sizeof(int)) ; x->d[0] = k+2 ; x->d[x->M] = k+2 ; /* sentinels */
}

void pbwtCursorCalculateU (PbwtCursor *x)
{
  int i, u = 0 ;
    
  for (i = 0 ; i < x->M ; ++i)
    { x->u[i] = u ;
      if (x->y[i] == 0)	++u ;
    }
  x->c = x->u[i] = u ;		/* need one off the end of update intervals */
}

/* We need to be careful about isBlockEnd in the routines below because when reading
   forwards we will have ->n naturally set after the end of the last block read,
   whereas when reading backwards it will naturally be at the start of the last block read.
   Also when writing it is at the start of the next block, because we don't yet have the contents.
*/

void pbwtCursorForwardsRead (PbwtCursor *u) /* move forwards and read (unless at end) */
{
  pbwtCursorForwardsAPacked (u) ;
  if (!u->isBlockEnd && u->n < arrayMax(u->z))  /* move to end of previous block */
    { u->nBlockStart = u->n ;
      u->n += unpack3 (arrp(u->z,u->n,uchar), u->M, u->y, 0) ;
    }
  if (u->n < arrayMax(u->z))
    { u->nBlockStart = u->n ;
      u->n += unpack3 (arrp(u->z,u->n,uchar), u->M, u->y, &u->c) ; /* read this block */
      u->isBlockEnd = TRUE ;
    }
  else
    u->isBlockEnd = FALSE ;	/* couldn't read in block and go to end */
}

void pbwtCursorForwardsReadAD (PbwtCursor *u, int k) /* AD version of the above */
{
  pbwtCursorForwardsAD (u, k) ;
  if (!u->isBlockEnd && u->n < arrayMax(u->z))  /* move to end of previous block */
    { u->nBlockStart = u->n ;
      u->n += unpack3 (arrp(u->z,u->n,uchar), u->M, u->y, 0) ;
    }
  if (u->n < arrayMax(u->z))
    { u->nBlockStart = u->n ;
      u->n += unpack3 (arrp(u->z,u->n,uchar), u->M, u->y, &u->c) ; /* read this block */
      u->isBlockEnd = TRUE ;
    }
  else
    u->isBlockEnd = FALSE ;	/* couldn't read in block and go to end */
}

void pbwtCursorReadBackwards (PbwtCursor *u) /* read and go backwards (unless at start) */
{
  if (u->isBlockEnd && u->n) u->n -= packCountReverse (arrp(u->z,u->n,uchar), u->M) ;
  if (u->n)
    { u->n -= packCountReverse (arrp(u->z,u->n,uchar), u->M) ;
      u->nBlockStart = u->n ;
      unpack3 (arrp(u->z,u->n,uchar), u->M, u->y, &u->c) ;
      pbwtCursorBackwardsA (u) ;
      u->isBlockEnd = FALSE ;
    }
  else
    u->isBlockEnd = TRUE ;
}

void pbwtCursorWriteForwards (PbwtCursor *u) /* write then move forwards */
{
  u->n += pack3arrayAdd (u->y, u->M, u->z) ;
  u->isBlockEnd = FALSE ;
  pbwtCursorForwardsA (u) ;
}

void pbwtCursorWriteForwardsAD (PbwtCursor *u, int k)
{
  u->n += pack3arrayAdd (u->y, u->M, u->z) ;
  u->isBlockEnd = FALSE ;
  pbwtCursorForwardsAD (u, k) ;
}

void pbwtCursorToAFend (PbwtCursor *u, PBWT *p) /* utility to copy final u->a to p->aFend */
{
  if (!p->aFend) p->aFend = myalloc (p->M, int) ; 
  memcpy (p->aFend, u->a, p->M*sizeof(int)) ;
}

/***************************************************/

void pbwtCursorForwardsAPacked (PbwtCursor *u)
/* A replacement for pbwtCursorForwardsA()
   We need u->nBlockStart = start of the packed array corresponding to current y,
   then copy blocks of old a into new a.
   Unfortunately we can't do this with AD because we need the maximum of d in a block.
*/
{
  int c = 0, m = 0, n ;
  uchar *zp = arrp(u->z,u->nBlockStart,uchar), *zp0 = zp, z ;
  while (m < u->M) {
    z = *zp++ ;
    n = p3decode[z & 0x7f] ; z >>= 7 ;
    if (z)
      memcpy (u->b+(m-c), u->a+m, n*sizeof(int)) ;
    else {
      memcpy (u->e+c, u->a+m, n*sizeof(int)) ;
      c += n ;
    }
    m += n ;
  }
  if (m != u->M) die ("error in forwardsAPacked()") ;

  memcpy (u->a, u->e, u->c*sizeof(int)) ;
  memcpy (u->a+u->c, u->b, (u->M-u->c)*sizeof(int)) ;
}


/***************************************************/

PBWT *pbwtSelectSites (PBWT *pOld, Array sites, BOOL isKeepOld)
{
  PBWT *pNew = pbwtCreate (pOld->M, 0) ;
  int ip = 0, ia = 0, j ;
  Site *sp = arrp(pOld->sites,ip,Site), *sa = arrp(sites,ia,Site) ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  uchar *x = myalloc (pNew->M, uchar) ;

  int nMissingSites = 0 ; 
  uchar *xMissing = myalloc (pNew->M+1, uchar) ;
  xMissing[pNew->M] = Y_SENTINEL ;  /* needed for efficient packing */

  pNew->sites = arrayCreate (arrayMax(sites), Site) ;
  while (ip < pOld->N && ia < arrayMax(sites))
    { if (sp->x < sa->x) 
        { ++ip ; ++sp ;
          pbwtCursorForwardsRead (uOld) ;
        }
      else if (sp->x > sa->x) { ++ia ; ++sa ; }
      else 
        {
          char *sa_als = dictName(variationDict, sa->varD) ;
          char *sp_als = dictName(variationDict, sp->varD) ;
          BOOL noAlt = sa_als[strlen(sa_als)-1] == '.' || sp_als[strlen(sp_als)-1] == '.';
          if (!noAlt && sp->varD < sa->varD)
            { ++ip ; ++sp ;
              pbwtCursorForwardsRead (uOld) ;
            }
          else if (!noAlt && sp->varD > sa->varD) { ++ia ; ++sa ; }
          else
            { array(pNew->sites,pNew->N,Site) = *sp ;
              for (j = 0 ; j < pOld->M ; ++j) x[uOld->a[j]] = uOld->y[j] ;
              pbwtCursorForwardsRead (uOld) ;
              for (j = 0 ; j < pNew->M ; ++j) uNew->y[j] = x[uNew->a[j]] ;
              pbwtCursorWriteForwards (uNew) ;

              if (pOld->missingOffset)
                { if (!pNew->missingOffset)
                    { pNew->zMissing = arrayCreate (10000, uchar) ;
                      array(pNew->zMissing, 0, uchar) = 0 ; /* needed so missing[] has offset > 0 */
                      pNew->missingOffset = arrayCreate (1024, long) ;
                    }
                  if (!arr(pOld->missingOffset,ip,long)) bzero (xMissing, pNew->M) ;
                  else unpack3 (arrp(pOld->zMissing, arr(pOld->missingOffset,ip,long), uchar), pNew->M, xMissing, 0) ;
                  if (arr(pOld->missingOffset,ip,long))
                    { 
                      array(pNew->missingOffset,pNew->N,long) = arrayMax(pNew->zMissing) ;
                      pack3arrayAdd (xMissing,pNew->M,pNew->zMissing) ; /* NB original order, not pbwt sort */
                      nMissingSites++ ;
                    }
                  else
                    array(pNew->missingOffset,pNew->N,long) = 0 ;
                }
              ++ip ; ++sp ; ++ia ; ++sa ; pNew->N++ ;
            }
        }
    }
  pbwtCursorToAFend (uNew, pNew) ;
  free(xMissing) ;

  // TODO: subset dosage

  fprintf (logFile, "%d sites selected from %d, %d missing sites, pbwt size for %d haplotypes is %ld\n", 
	   pNew->N, pOld->N, nMissingSites, pNew->M, arrayMax(pNew->yz)) ;

  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;

  if (isKeepOld)
    { if (pOld->samples) pNew->samples = arrayCopy (pOld->samples) ;
      if (pOld->chrom) pNew->chrom = strdup (pOld->chrom) ;
    }
  else				/* destroy one or the other */
    if (pNew->N == pOld->N)	/* no change - keep pOld as pNew */
      { pbwtDestroy (pNew) ;
	pNew = pOld ;
      }
    else
      { pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
	pNew->samples = pOld->samples ; pOld->samples = 0 ;
        if (pOld->missingOffset) arrayDestroy (pOld->missingOffset) ;
        pOld->missingOffset = 0 ;
        if (pOld->zMissing) arrayDestroy (pOld->zMissing) ;
        pOld->zMissing = 0 ;
	pbwtDestroy (pOld) ;
      }

  free(x) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  return pNew ;
}

/***************************************************/

PBWT *pbwtRemoveSites (PBWT *pOld, Array sites, BOOL isKeepOld)
{
  PBWT *pNew = pbwtCreate (pOld->M, 0) ;
  int ip = 0, ia = 0, j ;
  Site *sp = arrp(pOld->sites,ip,Site), *sa = arrp(sites,ia,Site) ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  uchar *x = myalloc (pNew->M, uchar) ;

  int nMissingSites = 0 ; 
  uchar *xMissing = myalloc (pNew->M+1, uchar) ;
  xMissing[pNew->M] = Y_SENTINEL ;  /* needed for efficient packing */

  pNew->sites = arrayCreate (arrayMax(sites), Site) ;
  while (ip < pOld->N)
    { if (ia>=arrayMax(sites) || sp->x < sa->x || (sp->varD < sa->varD && sp->x <= sa->x)) 
        { array(pNew->sites,pNew->N,Site) = *sp ;
          for (j = 0 ; j < pOld->M ; ++j) x[uOld->a[j]] = uOld->y[j] ;
          pbwtCursorForwardsRead (uOld) ;
          for (j = 0 ; j < pNew->M ; ++j) uNew->y[j] = x[uNew->a[j]] ;
          pbwtCursorWriteForwards (uNew) ;

          if (pOld->missingOffset)
            { if (!pNew->missingOffset)
                { pNew->zMissing = arrayCreate (10000, uchar) ;
                  array(pNew->zMissing, 0, uchar) = 0 ; /* needed so missing[] has offset > 0 */
                  pNew->missingOffset = arrayCreate (1024, long) ;
                }
              if (!arr(pOld->missingOffset,ip,long)) bzero (xMissing, pNew->M) ;
              else unpack3 (arrp(pOld->zMissing, arr(pOld->missingOffset,ip,long), uchar), pNew->M, xMissing, 0) ;
              if (arr(pOld->missingOffset,ip,long))
                { 
                  array(pNew->missingOffset,pNew->N,long) = arrayMax(pNew->zMissing) ;
                  pack3arrayAdd (xMissing,pNew->M,pNew->zMissing) ; /* NB original order, not pbwt sort */
                  nMissingSites++ ;
                }
              else
                array(pNew->missingOffset,pNew->N,long) = 0 ;
            }
          ++ip ; ++sp ; pNew->N++ ;
        }
      else if (sp->x > sa->x || sp->varD > sa->varD) { ++ia ; ++sa ; }
      else
	{ ++ip ; ++sp ; ++ia ; ++sa ;
	  pbwtCursorForwardsRead (uOld) ;
	}
    }
  pbwtCursorToAFend (uNew, pNew) ;

  fprintf (logFile, "%d sites selected from %d, %d missing sites, pbwt size for %d haplotypes is %ld\n", 
	   pNew->N, pOld->N, nMissingSites, pNew->M, arrayMax(pNew->yz)) ;

  if (isKeepOld)
    { if (pOld->samples) pNew->samples = arrayCopy (pOld->samples) ;
      if (pOld->chrom) pNew->chrom = strdup (pOld->chrom) ;
    }
  else				/* destroy one or the other */
    if (pNew->N == pOld->N)	/* no change - keep pOld as pNew */
      { pbwtDestroy (pNew) ;
	pNew = pOld ;
      }
    else
      { pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
	pNew->samples = pOld->samples ; pOld->samples = 0 ;
	pbwtDestroy (pOld) ;
      }
  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;

  free(x) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  return pNew ;
}

/******************* end of file *******************/

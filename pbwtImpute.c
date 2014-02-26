/*  File: pbwtImpute.c
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
 * Description: phasing and imputation functions in pbwt package, 
                plus utilities to intentionally corrupt data
 * Exported functions:
 * HISTORY:
 * Last edited: Feb  5 14:50 2014 (rd)
 * Created: Thu Apr  4 12:02:56 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

static void genotypeComparePbwt (PBWT *p, PBWT *q) ; /* forward declaration */

/************************* IMPUTATION AND PHASING *************************/

#include <math.h>

static double fBound[] = {0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100.01} ;

void imputeExplore (PBWT *p, int test)
{
  int i, j, k, M = p->M, N = p->N, ff ;
  uchar *y ;
  double xbar, ybar, tot, r2 ;
  double f ;
  typedef struct {
    long n00, n01 ;		/* neither neighbour is 1, truth is 0 or 1 */
    long n10a, n10b, n11a, n11b ;/* one neighbour is 1; a if d is lower for the 0 neighbour */
    long n20, n21 ;		/* both neighbours are 1 */
  } TestStat ;
  TestStat *testStat = mycalloc (16, TestStat), *t ;
  typedef long Counts[4] ;
  Array dHist = arrayCreate (1000, Counts) ;
  Counts cSimple, cCond0, cCond1 ;
  int *n0 = myalloc (M, int), *n1 = myalloc (M, int) ;
  uchar *x = myalloc (M, uchar) ;
  static long c0[16][5], c1[16][5] ;

  pbwtBuildReverse (p) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  PbwtCursor *uz = pbwtCursorCreate (p, FALSE, FALSE) ;

  for (i = 0 ; i < 4 ; ++i) cSimple[i] = cCond0[i] = cCond1[i] = 0 ; 

  for (k = 0 ; k < N ; ++k)
    { pbwtCursorReadBackwards (uz) ;
      if (isCheck)
	{ for (i = 0 ; i < M ; ++i)
	    x[u->a[i]] = u->y[i] ;
	  for (i = 0 ; i < M ; ++i)
	    if (x[uz->a[i]] != uz->y[i]) 
	      fprintf (stderr, "forward-backward mismatch at k %d i %d\n", k, i) ;
	}
      if (k > 0.2*N && k < 0.8*N)     /* ignore ends */
	{ f = (M - u->c) / (double)M ; for (ff = 0 ; f*100 > fBound[ff] ; ++ff) ;
	  t = &testStat[ff] ;
	  memset (n0, 0, M*sizeof(int)) ; memset (n1, 0, M*sizeof(int)) ; 
	  for (i = 1 ; i < M-1 ; ++i)
	    { if (u->y[i-1] && u->y[i+1])
		if (u->y[i]) ++t->n21 ; else ++t->n20 ;
	      else if (!u->y[i-1] && !u->y[i+1])
		if (u->y[i]) ++t->n01 ; else ++t->n00 ;
	      else if (!u->y[i-1] && u->d[i] < u->d[i+1] || !u->y[i+1] && u->d[i+1] < u->d[i])
		if (u->y[i]) ++t->n11a ; else ++t->n10a ;
	      else
		if (u->y[i]) ++t->n11b ; else ++t->n10b ;
	      ++(array(dHist, u->d[i]/100, Counts)[u->y[i-1] + 2*u->y[i]]) ;
	      ++cSimple[u->y[i-1] + 2*u->y[i]] ;
	      if (u->y[i+1])
		++cCond1[u->y[i-1] + 2*u->y[i]] ;
	      else
		++cCond0[u->y[i-1] + 2*u->y[i]] ;
	      n0[u->a[i]] += 2 - (u->y[i-1] + u->y[i+1]) ;
	      n1[u->a[i]] += u->y[i-1] + u->y[i+1] ;
	      n0[uz->a[i]] += 2 - (uz->y[i-1] + uz->y[i+1]) ;
	      n1[uz->a[i]] += uz->y[i-1] + uz->y[i+1] ;
	      x[u->a[i]] = u->y[i] ;
	    }
	  for (i = 0 ; i < M ; ++i)
	    if (n0[i] + n1[i] == 4)	/* it wasn't at the end either forwards or backwards */
	      if (x[i]) ++c1[ff][n1[i]] ; else ++c0[ff][n1[i]] ;
	}
      pbwtCursorForwardsReadADU (u, k) ;
    }

  if (test == 1)
    for (j = 0 ; j < 16 ; ++j)
      { t = &testStat[j] ;
	printf ("%-5.1f\t00,01\t%ld\t%ld\t10a,11a\t%ld\t%ld\t10b,11b\t%ld\t%ld\t20,21\t%ld\t%ld",
		fBound[j], t->n00, t->n01, t->n10a, t->n11a, t->n10b, t->n11b, t->n20, t->n21) ;
	tot = t->n00 + t->n01 + t->n10a + t->n11a + t->n10b + t->n11b + t->n20 + t->n21 ;
	if (tot)
	  { xbar = (t->n10b + t->n11b + t->n20 + t->n21)/tot ;
	    ybar = (t->n01 + t->n11a + t->n11b + t->n21)/tot ;
	    r2 = ((t->n21 + t->n11b)/tot - xbar*ybar)/sqrt((xbar - xbar*xbar)*(ybar - ybar*ybar)) ;
	    printf ("\tx,y,r2\t%.4f\t%.4f\t%.4f\n", xbar, ybar, r2) ;
	  }
	else
	  putchar ('\n') ;
      }
  else if (test == 2)
    for (j = 0 ; j < arrayMax(dHist) ; ++j)
      { long *c ; c = arr(dHist, j, Counts) ;
	printf ("%d\t%ld\t%ld\t%ld\t%ld", j, c[0], c[1], c[2], c[3]) ;
	if (c[0] + c[2]) printf ("\t%.3f", c[0]/(double)(c[0]+c[2])) ; else printf ("\t0") ;
	if (c[1] + c[3]) printf ("\t%.3f", c[3]/(double)(c[1]+c[3])) ; else printf ("\t0") ;
	putchar ('\n') ;
      }
  else if (test == 3)
    { printf ("%.3f %.3f\t", cSimple[0]/(double)(cSimple[0]+cSimple[2]), cSimple[3]/(double)(cSimple[1]+cSimple[3])) ;
      printf ("%.3f %.3f\t", cCond0[0]/(double)(cCond0[0]+cCond0[2]), cCond0[3]/(double)(cCond0[1]+cCond0[3])) ;
      printf ("%.3f %.3f\n", cCond1[0]/(double)(cCond1[0]+cCond1[2]), cCond1[3]/(double)(cCond1[1]+cCond1[3])) ;
    }
  else if (test == 4)
    for (j = 0 ; j < 16 ; ++j)
      { printf ("%-5.1f", fBound[j]) ;
	tot = 0 ; xbar = 0 ; r2 = 0 ;
        for (i = 0 ; i < 5 ; ++i)
	  { long sum = c0[j][i] + c1[j][i] ;
	    printf ("\t%ld ", sum) ;
	    if (sum) printf (" %.3f", c1[j][i]/(double)sum) ; else printf (" 00000") ;
	    xbar += c1[j][i] ;
	    if (i == 3 || i == 4) r2 += c1[j][i] ; if (i ==2) r2 += 0.5*c1[j][i] ;
	    tot += sum ;
	  }
	ybar = c0[j][4] + c1[j][4] + c0[j][3] + c1[j][3] + 0.5*(c0[j][2] + c1[j][2]) ;
	if (tot)
	  { xbar /= tot ;
	    ybar /= tot ;
	    r2 = (r2/tot - xbar*ybar)/sqrt((xbar - xbar*xbar)*(ybar - ybar*ybar)) ;
	    printf ("\tx,y,r2\t%.4f\t%.4f\t%.4f\n", xbar, ybar, r2) ;
	  }
	else
	  putchar ('\n') ;
      }

  pbwtCursorDestroy (u) ;
  pbwtCursorDestroy (uz) ;
}

/****************** phasing ****************/

static void phaseCompare (PBWT *p, PBWT *q)
{
  int i, k ;
  int M = p->M, N = p->N ;
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  PbwtCursor *uq = pbwtCursorCreate (q, TRUE, TRUE) ;
  int *xp = myalloc (M, int), *xq = myalloc (M, int) ;
  int *isFirst = myalloc (M, int), *isFlipped = myalloc (M, int) ;
  int *lastFlip = mycalloc (M, int), *kHet = mycalloc (M, int) ;
  int nSwitch = 0, nHet = 0, nSwitch1 = 0, nSwitch5 = 0 ;
  double mFac = 2.0/M ;
  int *nSwitchSample = mycalloc (M, int), *nSwitchSite = mycalloc (N, int) ;

  if (p->M != q->M || p->N != q->N) die ("size incompatibility in phaseCompare") ;
  if (p->M %2) die ("phaseCompare requires that M %d is even", M) ;

  memset (isFirst, 1, M*sizeof(int)) ; /* I know this is not 1 - anything non-zero will do */
  for (k = 0 ; k < N ; k++)
    { for (i = 0 ; i < M ; ++i)
	{ xp[up->a[i]] = up->y[i] ;
	  xq[uq->a[i]] = uq->y[i] ;
	}
      for (i = 0 ; i < M ; i += 2)
	{ if (xp[i] + xp[i+1] == 1)
	    { ++nHet ; ++kHet[i] ;
	      if (isFirst[i])
		{ isFirst[i] = 0 ;
		  isFlipped[i] = (xp[i] == xq[i+1]) ? 1 : 0 ; 
		}
	      else if (xp[i] != xq[i+isFlipped[i]]) 
		{ ++nSwitch ; ++nSwitchSample[i/2] ; ++nSwitchSite[k] ;
		  if (kHet[i] - lastFlip[i] > 1) ++nSwitch1 ;
		  if (kHet[i] - lastFlip[i] > 5) ++nSwitch5 ;
		  isFlipped[i] = 1 - isFlipped[i] ;
		  lastFlip[i] = kHet[i] ;
		}
	    }
	  if (isCheck && (xp[i]+xp[i+1] != xq[i]+xq[i+1])) die ("phaseCompare mismatch: k %d, i %d", k, i) ;
	}
      pbwtCursorForwardsRead (up) ;
      pbwtCursorForwardsRead (uq) ;
    }

  fprintf (stderr, "%.1f switches per sample, %.3f per het, %.1f nSwitch1, %.1f nSwitch5\n", 
	   mFac*nSwitch, nSwitch/(double)nHet, mFac*nSwitch1, mFac*nSwitch5) ;

  if (isStats)
    { for (i = 0 ; i < M/2 ; ++i)
	{ printf ("SAMPLE-SWITCH\t%d\t%d", i, nSwitchSample[i]) ;
	  if (p->samples)
	    printf ("\t%s", sampleName(arr(p->samples, 2*i, int))) ;
	  putchar ('\n') ;
	}
      for (k = 0 ; k < N ; ++k)	
	{ printf ("SITE-SWITCH\t%d\t%d", k, nSwitchSite[k]) ;
	  if (p->sites)
	    { Site *s = arrp(p->sites,k,Site) ;
	      printf ("\t%s\t%d\t%s", p->chrom, s->x, dictName(variationDict, s->varD)) ;
	    }
	  putchar ('\n') ;
	}
    }

  pbwtCursorDestroy (up) ; pbwtCursorDestroy (uq) ;
  free (xp) ; free (xq) ; free (isFirst) ; free (isFlipped) ; 
  free (nSwitchSample) ;  free (nSwitchSite) ;
}

/****************************************************************/

static double *scoreBit ;
static double *logisticCache ;

static void phaseInit (int N)
{ int i ;
  double z ;
  scoreBit = myalloc (N+1, double) ;
  for (i = 0 ; i <= N ; ++i) scoreBit[i] = log (i + 1.0) ; /* 1 is simple version */
  logisticCache = myalloc (100000, double) ;
  for (i = 0 ; i < 100000 ; ++i) 
    { z = exp (-i * 0.0001) ; logisticCache[i] = 1.0 / (1.0 + z) ; }
}

static inline double score0 (PbwtCursor *u, double *xp, int i)
{
  double s = 0.0 ;
  int ubi = u->b[i] ;
  if (ubi > 0) s += xp[u->a[ubi-1]] ;
  if (ubi < u->M-1) s += xp[u->a[ubi+1]] ;
  return s ;
}

static inline double score1 (PbwtCursor *u, double *xp, int i, int k)
{
  double s = 0 ;
  int ubi = u->b[i] ;
  if (ubi > 0) s += xp[u->a[ubi-1]] * scoreBit[(k+1)-u->d[ubi]] ;
  if (ubi < u->M-1) s += xp[u->a[ubi+1]] * scoreBit[(k+1)-u->d[ubi+1]] ;
  return s ;
}

static inline double logistic (double x)
{ int i = x * 10000 ; 
  if (i < -99999) return 1.0 - logisticCache[99999] ;
  else if (i < 0) return 1.0 - logisticCache[-i] ;
  else if (i < 100000) return logisticCache[i] ;
  else return logisticCache[99999] ;
}

static int initialiseXpFromX (double *xp, uchar *x, int M)
{
  int i, n2 = 0 ;
  for (i = 0 ; i < M ; i += 2) /* initialise xp */
    if (x[i] != x[i+1]) /* a het */
      { xp[i] = xp[i+1] = 0.0 ; ++n2 ; }
    else if (x[i] == 1) /* && x[i+1] == 1 */
      { xp[i] = xp[i+1] = 1.0 ; }
    else /* x[i] == x[i+1] == 0 */
      { xp[i] = xp[i+1] = -1.0 ; }
  return n2 ;
}

PBWT *phase (PBWT *p, int kMethod, int nSparse) /* return rephased p */
{
  int    i, j, k, kk, kr = 0 ;
  int    M = p->M, N = p->N ;
  PbwtCursor *up = pbwtNakedCursorCreate (M, 0) ;
  PBWT   *q ;		            /* new pbwt */
  PbwtCursor *uq, *ur, **uqq, **urr ;   /* update objects for forward and reverse, uqq and urr for sparse */
  int    nyp = 0, nyr, ny, n2, n2Old ;
  uchar  *x = myalloc (M, uchar) ;  /* actual haplotypes in original order, from p */
  double *xp = myalloc(M, double) ; /* 2*p(x=1)-1, so 1 if x=1, -1 if x=0, 0 if unknown */
  BOOL   is2way = (kMethod > 2) ? TRUE : FALSE ;
  BOOL   *isReverseSite = myalloc (N, BOOL) ;
  
  if (!p) die ("no pbwt for phase") ;
  if (p->M %2) die ("phaseTest requires that M = %d is even", M) ;
  if (kMethod < 0 || kMethod > 4 || nSparse < 0)
    die ("kMethod %d or nSparse %d out of range in phase", kMethod, nSparse) ;

  if (nSparse == 1) nSparse = 0 ; /* no point in just one nSparse */

  phaseInit (N) ;		/* sets up some lookup tables */

  q = pbwtCreate (M, N) ;

  if (is2way)		/* build reverse model first */
    { double goodSiteThresh = 1.0 ;
      int nSitesUsed = 0 ;

      /* first initialise reverse pbwt and associated data structures */
      q->zz = arrayCreate (arrayMax(q->yz), uchar) ;
      ur = pbwtNakedCursorCreate (M, 0) ;
      for (i = 0 ; i < M ; ++i) ur->b[ur->a[i]] = i ;  /* use u->b for inverse of u->a */
      if (nSparse)
	{ urr = myalloc (nSparse, PbwtCursor*) ;
	  for (kk = 0 ; kk < nSparse ; ++kk)
	    { urr[kk] = pbwtNakedCursorCreate (M, 0) ; 
	      for (i = 0 ; i < M ; ++i) urr[kk]->b[urr[kk]->a[i]] = i ;
	    }
	}

      for (i = 0 ; i < p->N ; ++i)	/* first run forwards to the end to set up */
	{ nyp += unpack3 (arrp(p->yz,nyp,uchar), M, up->y, 0) ;
	  pbwtCursorForwardsA (up) ;
	}
      fprintf (stderr, "Forward 1: ") ; timeUpdate () ;

      /* run back through p, building ur */
      for (k = N ; k-- ; )
	{ nyp -= packCountReverse (arrp(p->yz, nyp, uchar), M) ;
	  unpack3 (arrp(p->yz, nyp, uchar), M, up->y, &up->c) ;
	  pbwtCursorBackwardsA (up) ; /* this needs to come before the next section */
	  for (i = 0 ; i < M ; ++i) x[up->a[i]] = up->y[i] ;  /* build x from up->y */
	  n2 = initialiseXpFromX (xp, x, M) ;

	  { double thresh = 2*(nSparse+1) + 0.5 ;
	    double s ;

	    n2Old = M+1 ;
	    while (n2 && n2 < n2Old)   /* stop if no longer reducing unresolved genotypes */
	      { n2Old = n2 ; n2 = 0 ;
		for (i = 0 ; i < M ; i += 2) /* loop over genotype pairs in original order */
		  if (!xp[i]) /* a het to phase */
		    { s = score0 (ur, xp, i) - score0 (ur, xp, i+1) ;
		      for (kk = 0 ; kk < nSparse ; ++kk) 
			s += score0 (urr[kk], xp, i) - score0 (urr[kk], xp, i+1) ;
		      if (s > thresh)  { xp[i] = 1 ; xp[i+1] = -1 ; }
		      else if (s < -thresh) { xp[i] = -1 ; xp[i+1] = 1 ; }
		      else ++n2 ;
		    }
		if (n2 == n2Old && thresh > 1.0) { thresh -= 1.0 ; n2Old = M+1 ; }
	      }
	    if (n2)   /* some unresolved values - phase using length */
	      for (i = 0 ; i < M ; i += 2)
	      if (!xp[i]) 
		/* a little care needed here: we need to use reverse coords for ur */
		{ s = score1 (ur, xp, i, N-k) - score1 (ur, xp, i+1, N-k) ;
		  for (kk = 0 ; kk < nSparse ; ++kk) 
		    s += score1 (urr[kk], xp, i, (N-k)/nSparse) - score1 (urr[kk], xp, i+1, (N-k)/nSparse) ;
		  if (s > 0) { xp[i] = 1 ; xp[i+1] = -1 ; }
		  else { xp[i] = -1 ; xp[i+1] = 1 ; }
		}
	  }

	  /* choose final phasing */
	  for (i = 0 ; i < M ; ++i) x[i] = (xp[i] > 0.0) ? 1 : 0 ;

	  if (kMethod == 4)	/* decide whether to use this site or not */
	    { int nFlip = 0 ; double flipRate ;
	      int c = up->c ;
	      if (2*c > M) c = M - c ;
	      for (i = 1 ; i < M ; ++i) if (x[ur->a[i]] != x[ur->a[i-1]]) ++nFlip ;
	      for (kk = 0 ; kk < nSparse ; ++kk)
		for (i = 1 ; i < M ; ++i) if (x[urr[kk]->a[i]] != x[urr[kk]->a[i-1]]) ++nFlip ;
	      flipRate = nFlip/(double)((1+nSparse)*2*c) ;
	      isReverseSite[k] = (flipRate < goodSiteThresh) ;
	      /* printf ("nFlip %d, c %d, flipRate %.3f, goodSiteThresh %.2f\n",
		 nFlip, c, flipRate, goodSiteThresh) ; */
	      if (isReverseSite[k])
		goodSiteThresh -= 0.01 ;
	      else
		goodSiteThresh += 0.03 ;
	    }
	  else
	    isReverseSite[k] = TRUE ;

	  if (isReverseSite[k])
	    { ++nSitesUsed ;

	      /* update reverse pbwt */
	      for (i = 0 ; i < M ; ++i) ur->y[i] = x[ur->a[i]] ;
	      pack3arrayAdd (ur->y, M, q->zz) ;

	      /* and necessary update structures */
	      pbwtCursorForwardsADU (ur, nSitesUsed) ; 
	      for (i = 0 ; i < M ; ++i) ur->b[ur->a[i]] = i ;
	      if (nSparse)
		{ kk = k % nSparse ;	/* which of the sparse pbwts to update this time */
		  for (i = 0 ; i < M ; ++i) urr[kk]->y[i] = x[urr[kk]->a[i]] ;
		  pbwtCursorForwardsADU (urr[kk], (N-k)/nSparse) ; 
		  for (i = 0 ; i < M ; ++i) urr[kk]->b[urr[kk]->a[i]] = i ;
		}
	    }
	}


      if (isCheck)		/* flip uq->zz round into r and compare to p */
	{ PBWT *r = pbwtCreate (M, N) ;
	  r->yz = q->zz ; pbwtBuildReverse (r) ; r->yz = r->zz ; r->zz = 0 ;
	  fprintf (stderr, "After reverse pass: ") ; phaseCompare (p, r) ;
	  pbwtDestroy (r) ;
	}

      fprintf (stderr, "Reverse pass: %d of %d sites used, final thresh %.2f\n", 
	       nSitesUsed, N, goodSiteThresh) ;
      timeUpdate () ;
    }

  /* now build forward data structures */
  uq = pbwtNakedCursorCreate (M, 0) ; 
  /* if (is2way) memcpy (uq->a, ur->a, M*sizeof(int)) ; /* prime uq with final ur */
  /* why was the previous line wrong - it clearly did not work */
  for (i = 0 ; i < M ; ++i) uq->b[uq->a[i]] = i ;
  if (nSparse)
    { uqq = myalloc (nSparse, PbwtCursor*) ;
      for (kk = 0 ; kk < nSparse ; ++kk)
	{ uqq[kk] = pbwtNakedCursorCreate (M, 0) ; 
	  if (is2way) memcpy (uqq[kk]->a, urr[kk]->a, M*sizeof(int)) ;
	  for (i = 0 ; i < M ; ++i) uqq[kk]->b[uqq[kk]->a[i]] = i ;
	}
    }

  /* and run forwards to build the final phasing */
  nyp = 0 ; if (is2way) nyr = arrayMax(q->zz) ;
  for (k = 0 ; k < N ; k++)
    { nyp += unpack3 (arrp(p->yz,nyp,uchar), M, up->y, &up->c) ;
      for (i = 0 ; i < M ; ++i) x[up->a[i]] = up->y[i] ;  /* build x from up->y */
      pbwtCursorForwardsA (up) ;
      n2 = initialiseXpFromX (xp, x, M) ;

      /* method 0: iterative majority voting on neighbours - fast and surprisingly good */
      if (kMethod == 0 || is2way)
	{ double s ;
	  double thresh = 2*(nSparse + is2way?2:1) + 0.5 ;

	  n2Old = M+1 ;
	  while (n2 && n2 < n2Old)   /* stop if no longer reducing unresolved genotypes */
	    { n2Old = n2 ; n2 = 0 ;
	      for (i = 0 ; i < M ; i += 2) /* loop over genotype pairs in original order */
		if (!xp[i]) /* a het to phase */
		  { s = score0 (uq, xp, i) - score0 (uq, xp, i+1) ;
		    if (is2way) s += score0 (ur, xp,i) - score0 (ur, xp, i+1) ;
		    for (kk = 0 ; kk < nSparse ; ++kk) 
		      { s += score0 (uqq[kk], xp, i) - score0 (uqq[kk], xp, i+1) ;
#ifdef NO_HELP
			if (is2way) s += score0 (urr[kk], xp, i) - score0 (urr[kk], xp, i+1) ;
#endif
		      }
		    if (s > thresh)  { xp[i] = 1 ; xp[i+1] = -1 ; }
		    else if (s < -thresh) { xp[i] = -1 ; xp[i+1] = 1 ; }
		    else ++n2 ;
		  }
	      if (n2 == n2Old && thresh > 1.0) { thresh -= 1.0 ; n2Old = M+1 ; }
	    }
	  if (n2)   /* some unresolved values - phase using length, forwards only for now */
	    for (i = 0 ; i < M ; i += 2)
	      if (!xp[i]) 
		{ s = score1 (uq, xp, i, k) - score1 (uq, xp, i+1, k) ;
		  for (kk = 0 ; kk < nSparse ; ++kk) 
		    s += score1 (uqq[kk], xp, i, k/nSparse) - score1 (uqq[kk], xp, i+1, k/nSparse) ;
		  if (s > 0) { xp[i] = 1 ; xp[i+1] = -1 ; }
		  else { xp[i] = -1 ; xp[i+1] = 1 ; }
		}
	}

      /* methods 1,2: iterative use of logistic sum of log-weighted shared prefix lengths  */
      else if (kMethod == 1 || kMethod == 2)	/* use logistic predictor */
	/* method 1 is mean field relaxation; method 2 is sampling - 2 a bit dissappointing */
	{ double lambda = 0.5 / (1+nSparse) ;	/* I tried fitting it, but fixed is better */

	  for (j = 0 ; j < 20 ; ++j) /* arbitrary limit of 20 iterations */
	    { n2 = 0 ;
	      for (i = 0 ; i < M ; i += 2) /* update xp[] for sites to phase */
		if (x[i] != x[i+1])	   /* need to phase */
		  { double s, newXPi ;
		    s = score1 (uq, xp, i, k) - score1 (uq, xp, i+1, k) ;
		    for (kk = 0 ; kk < nSparse ; ++kk) 
		      s += score1 (uqq[kk], xp, i, k/nSparse) - score1 (uqq[kk], xp, i+1, k/nSparse) ;
		    newXPi= 2 * logistic (lambda*s) - 1.0 ;
		    if (newXPi * xp[i] < 0.01) ++n2 ; /* xp[i] not set or sign changed */
		    if (kMethod == 1) /* expected value */
		      xp[i] = newXPi ; 
		    else  /* KMethod == 2: sample from posterior */
		      xp[i] = (rand() < 0.5*(newXPi+1)*RAND_MAX) ? 1.0 : -1.0 ;
		    xp[i+1] = -xp[i] ;
		  }
	      if (!n2) break ;
	    }
	}

      /* choose final phasing */
      for (i = 0 ; i < M ; ++i)	x[i] = (xp[i] > 0.0) ? 1 : 0 ;

      /* update q pbwt */
      for (i = 0 ; i < M ; ++i)	uq->y[i] = x[uq->a[i]] ;
      pack3arrayAdd (uq->y, M, q->yz) ;

      /* and related update structures */
      pbwtCursorForwardsADU (uq, k) ; for (i = 0 ; i < M ; ++i) uq->b[uq->a[i]] = i ;
      if (nSparse)
	{ kk = k % nSparse ;	/* which of the sparse pbwts to update this time */
	  for (i = 0 ; i < M ; ++i) uqq[kk]->y[i] = x[uqq[kk]->a[i]] ;
	  pbwtCursorForwardsADU (uqq[kk], k/nSparse) ; for (i = 0 ; i < M ; ++i) uqq[kk]->b[uqq[kk]->a[i]] = i ;
	}

      if (isStats)		/* report on this site */
	{ int c = up->c ;
	  if (2*c > M) c = M - c ;
	  printf ("SITE-INFO\t%d\tc\t%d", k, c) ;
	  if (c)
	    { int nFlip = 0 ; double fac = 0.5 / c ;
	      for (i = 1 ; i < M ; ++i) if (up->y[i] != up->y[i-1]) ++nFlip ;
	      printf ("\tp %.2f", nFlip*fac) ; nFlip = 0 ;
	      for (i = 1 ; i < M ; ++i) if (x[uq->a[i]] != x[uq->a[i-1]]) ++nFlip ;
	      printf ("\tq %.2f", nFlip*fac) ; nFlip = 0 ;
	      for (kk = 0 ; kk < nSparse ; ++kk)
		{ for (i = 1 ; i < M ; ++i) 
		    if (x[uqq[kk]->a[i]] != x[uqq[kk]->a[i-1]]) ++nFlip ;
		  printf ("\tq%d %.2f", kk, nFlip*fac) ; nFlip = 0 ;
		}
	      if (is2way)
		{ for (i = 1 ; i < M ; ++i) if (x[ur->a[i]] != x[ur->a[i-1]]) ++nFlip ;
		  printf ("\tr %.2f", nFlip*fac) ; nFlip = 0 ;
#ifdef NO_HELP
		  for (kk = 0 ; kk < nSparse ; ++kk)
		    { nFlipR = 0 ;
		      for (i = 1 ; i < M ; ++i) 
			if (x[urr[kk]->a[i]] != x[urr[kk]->a[i-1]]) ++nFlipR ;
		      printf ("\t%.3f", 0.5*nFlipR/(double)c) ;
		    }
#endif
		}
	    }
	  putchar ('\n') ;
	}

      /* if we have a reverse pbwt, move ur and urr back in that */
      if (is2way && isReverseSite[k])
	{ nyr -= packCountReverse (arrp(q->zz, nyr, uchar), M) ;
	  unpack3 (arrp(q->zz, nyr, uchar), M, ur->y, &ur->c) ;
	  for (i = 0 ; i < M ; ++i) x[ur->a[i]] = ur->y[i] ;  /* need x from ur, not uq, to undo urr[] */
	  pbwtCursorBackwardsA (ur) ; for (i = 0 ; i < M ; ++i) ur->b[ur->a[i]] = i ;
#ifdef NO_HELP
	  if (nSparse)
	    { kk = k % nSparse ;	/* which of the sparse pbwts to update this time */
	      for (i = 0 ; i < M ; ++i) urr[kk]->y[i] = x[urr[kk]->a[i]] ;
	      urr[kk]->c = ur->c ; 
	      pbwtCursorBackwardsA (urr[kk]) ;
	      for (i = 0 ; i < M ; ++i) urr[kk]->b[urr[kk]->a[i]] = i ;
	    }
#endif
	}
    }
  pbwtCursorToAFend (uq, q) ;

  /* destroy reverse PBWT, which is not the reverse of the forwards one at this point */
  free (q->zz) ; q->zz = 0 ;

  /* compare new phasing to original and report switch rates */
  fprintf (stderr, "After forward pass: ") ; phaseCompare (p, q) ;

  /* clean up memory allocated */
  free (x) ; free (xp) ;
  pbwtCursorDestroy (up) ; pbwtCursorDestroy (uq) ;
  if (nSparse) { for (kk = 0 ; kk < nSparse ; ++kk) free (uqq[kk]) ; free (uqq) ; }
  if (is2way) 
    { pbwtCursorDestroy (ur) ;
      if (nSparse) { for (kk = 0 ; kk < nSparse ; ++kk) free (urr[kk]) ; free (urr) ; }
    }
  pbwtDestroy (p) ;
  return q ;
}

/******* phase a new pbwt against the existing one as a reference *******/

typedef struct { double alpha, beta ; } ScoreParams ;

static inline ScoreParams scoreParams (PbwtCursor *u)
{
  ScoreParams sp ;
  sp.alpha = 0.0 ;		/* looks like for gen10k 100 is marginally better */
  sp.beta = 1.0 ;
  return sp ;
}

static ScoreParams scoreParamsTheory (PbwtCursor *u)
/* Vladimir says that p(0|1) = exp (-alpha - beta*d)  */
/* log(p) is about -alpha - beta*d
   log(1-p) is about -exp(-alpha - beta*d)
   let x_j = 1 if a change, else 0
   total log likelihood L = - sum_j x_j (alpha + beta*d_j) + (1-x_j) exp (-alpha - beta*d_j)
   dL/dalpha = 0 gives sum_j x_j = sum_j (1-x_j) exp (-alpha - beta*d_j)
   so alpha  = ln ((sum_j (1-x_j) exp (-beta*d_j))/sum_j x_j)
   dL/dbeta  = 0 gives sum_j x_j d_j = sum_j (1-x_j) d_j exp (-alpha - beta*d_j)
   substituting in alpha we get 
   sum_j x_j d_j / sum_j x_j = sum_j (1-x_j) d_j exp(-beta*d_j) / sum_j (1-x_j) exp(-beta*d_j)
*/
{
  ScoreParams sp ;
  int j ;
  double sum_dx = 0, sum_e, sum_de ; int sum_x = 0 ;

  if (u->c == 0 || u->c == u->M) { sp.alpha = 0 ; sp.beta = 0 ; return sp ; }

  for (j = 1 ; j < u->M ; ++j)
    if (u->y[j] != u->y[j-1])
      { sum_dx += u->d[j] ; ++sum_x ; }
  double T = sum_dx / sum_x ;	/* target */
  for (j = 1 ; j < u->M ; ++j)
    if (u->y[j] == u->y[j-1])
      ;

  sp.alpha = log (sum_e / sum_x) ;
      
  return sp ;
}

/************ routines to manage matches of sequences into pbwts ************/

typedef struct {
  int i ;    /* the position in the reference uRef->y that new sequence PRECEDES */
  int dminus, dplus ;	/* the start of the match to current ref sequences i-1 and i */
} MatchInfo ;

/* NB when allocating MatchInfo use mycalloc() so dminus, dplus are initialised to 0 */

static inline double matchPhaseScore (MatchInfo *m, PbwtCursor *u, ScoreParams *sp, int k)
{
  double score = 0 ;
  if (m->i > 0) score += (2*u->y[m->i-1] - 1) * (sp->alpha + sp->beta * (k - m->dminus)) ;
  if (m->i < u->M) score += (2*u->y[m->i] - 1) * (sp->alpha + sp->beta * (k - m->dplus)) ;
  ++m ;				/* move on to the second allele */
  if (m->i > 0) score -= (2*u->y[m->i-1] - 1) * (sp->alpha + sp->beta * (k - m->dminus)) ;
  if (m->i < u->M) score -= (2*u->y[m->i] - 1) * (sp->alpha + sp->beta * (k - m->dplus)) ;
  return score ;
}

static void matchUpdate (MatchInfo *match, uchar *xx, int M, PbwtCursor *u, int k)
/* This is ugly because it incorporates the pbwt cursor update on u in the middle of 
   updating match.  This is necessary because we need the old d[] but the new u[].
   We could make things tighter by integrating the match and cursor updates, perhaps
   even operating on the compressed space, but that is probably only viable for
   single rather than multiple matches.
*/
{
  int i, j ; MatchInfo *m ; uchar *x ;

  for (j = 0, m = match, x = xx ; j < M ; ++j, ++m, ++x)
    {  /* If x == y[i] then dplus does not change and we are done.
	  Otherwise, need to find next matching symbol.  Same in minus direction.
       */
      for (i = m->i ; i < u->M ; ++i)
	{ if (u->y[i] == *x) break ;
	  if (u->d[i] > m->dplus) m->dplus = u->d[i] ;
	}
      for (i = m->i-1 ; i >= 0 ; --i)
	{ if (u->y[i] == *x) break ;
	  if (u->d[i] > m->dminus) m->dminus = u->d[i] ;
	}
    }
  int c = u->c ;     		     /* need to cache this before updating u */
  pbwtCursorForwardsReadADU (u, k) ; /* need new u->u values to update m->i */
  for (j = 0, m = match, x = xx ; j < M ; ++j, ++m, ++x)
    m->i = *x ? c + m->i - u->u[m->i] : u->u[m->i] ; /* simple, isn't it! */

  if (isCheck)
    for (j = 0, m = match, x = xx ; j < M ; ++j, ++m, ++x)
       if (m->i < 0 || m->i > u->M) die ("out of bounds in matchUpdate") ;
}
/* Consistency checks on m->i update:
   m->i between 0 and M inclusive is the position just AFTER the query.
   u->u[i] is the number of 0s in y[] before position i.  
   So if x == 0 and there are no 0s before i = m->i then m->i is mapped to 0 - good.
   If x == 1 and there are no 1s then u[m->i] = m->i so m->i is mapped to c - good.
   If x == 1 and all 1s are before m->i then m->i is mapped to M - good.
   Mapping of 0 is correct, to 0 if 0 else c if 1 - good.
   Mapping of M is correct, to c if 0 else M if 1, so long as u[M] is defined as c - good.
   So it looks like this is all good.
*/

/************** main function to phase against a reference ***************/

static PBWT *referencePhase1 (PBWT *pOld, PBWT *pRef)
{
  if (!pRef->zz) pbwtBuildReverse (pRef) ;	/* we need the reverse reference pbwt below */

/* Now phase the new sites against the old sites
   We will pass forwards, then back, then forwards again.
   All we store on the first two passes is the phasing score at each ambiguous genotype from
   the left then right respectively, favouring 1/0 if the score is positive, and 0/1 if the score 
   is negative. 
*/
  /* declarations and initialisations */
  int M = pOld->M ;
  int j, k ;			    /* indices: new sample, site */
  PbwtCursor *uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;   /* cursor on old */
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;   /* cursor on new */
  uchar *x = myalloc (M, uchar) ;   /* current uOld values in original (unsorted) order */
  MatchInfo *match = mycalloc (M, MatchInfo) ;
  for (j = 0 ; j < M ; ++j)	    /* initialise m->i[] randomly */
    match[j].i = rand() * (double)pRef->M / RAND_MAX ;
  Array *qScoreStack = myalloc (M, Array) ; /* score from left at phase-ambiguous genotypes */
  for (j = 0 ; j < M ; j += 2)	    /* only need to store score for even positions */
    qScoreStack[j] = arrayCreate (1024, double) ; 
  ScoreParams sp ;

  /* first forwards loop */
  for (k = 0 ; k < pRef->N ; ++k)
    { for (j = 0 ; j < M ; ++j) /* extract genotype */
	x[uOld->a[j]] = uOld->y[j] ;
      sp = scoreParams (uRef) ;
      for (j = 0 ; j < M ; j += 2) /* step through in pairs and phase if a heterozygote */
	if (x[j] + x[j+1] == 1)
	  { double score = matchPhaseScore (&match[j], uRef, &sp, k) ;
	    /* for (kk = 0 ; kk < nSparse ; ++kk) score += matchPhaseScore (&matchSparse[kk][j], uRefSparse[kk], &sp, k/nSparse) ; */
	    array(qScoreStack[j], arrayMax(qScoreStack[j]), double) = score ; /* store here */
	    /* set phasing based on score - default when score == 0 to lexicographic order */
	    if (score > 0) { x[j] = 1 ; x[j+1] = 0 ; } else { x[j] = 0 ; x[j+1] = 1 ; }
	  }

      /* Update mapping into Ref: must do in two steps. Do this in a subroutine because
         we will do it several times.  NB the major side effect to move forwards in uRef */
      matchUpdate (match, x, M, uRef, k) ;
      /* if (nSparse) { kk = k % nsparse ; matchUpdate (matchSparse[kk], x, M, uRefSparse[k], k/nSparse) ; } */
      pbwtCursorForwardsRead (uOld) ;     /* finally move forwards in pOld */
    }
  fprintf (stderr, "first forward pass complete\n") ;

  /* Now reverse, carrying out the equivalent actions, but phasing from both directions.
     This time store the phasing score from the right.
  */
  Array *rScoreStack = myalloc (M, Array) ; /* stack for reverse scores at phase sites */
  for (j = 0 ; j < M ; j += 2)
    rScoreStack[j] = arrayCreate (arrayMax(qScoreStack[j]), double) ; 
  for (j = 0 ; j < M ; ++j) match[j].dminus = match[j].dplus = 0 ; /* reset d* */
  pbwtCursorDestroy (uRef) ; uRef = pbwtCursorCreate (pRef, FALSE, TRUE) ;
  /* note run uRef in the reverse pbwt - second arg of pbwtCursorCreate() is FALSE */

  for (k = 0 ; k < pRef->N ; ++k)
    { pbwtCursorReadBackwards (uOld) ;
      for (j = 0 ; j < M ; ++j) /* extract genotype */
	x[uOld->a[j]] = uOld->y[j] ;
      for (j = 0 ; j < M ; j += 2) /* step through in pairs and phase if necessary */
	if (x[j] + x[j+1] == 1)	 /* need to phase - add new right score to stored left score */
	  { double score = matchPhaseScore (&match[j], uRef, &sp, k) ;
	    array(rScoreStack[j], arrayMax(rScoreStack[j]), double) = score ; /* push onto rStoreStack */
	    score += arr(qScoreStack[j], --arrayMax(qScoreStack[j]), double) ; /* pop off qStoreStack */
	    if (score > 0) { x[j] = 1 ; x[j+1] = 0 ; } else { x[j] = 0 ; x[j+1] = 1 ; }
	  }
      matchUpdate (match, x, M, uRef, k) ;
    }
  fprintf (stderr, "reverse pass complete\n") ;

 /* Now the final forward pass, building pNew */
  for (j = 0 ; j < M ; ++j) match[j].dminus = match[j].dplus = 0 ; /* reset d* */
  pbwtCursorDestroy (uRef) ; uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  PBWT *pNew = pbwtCreate (M, pOld->N) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;

  for (k = 0 ; k < pRef->N ; ++k)
    { for (j = 0 ; j < M ; ++j) /* extract genotype */
	x[uOld->a[j]] = uOld->y[j] ;
      for (j = 0 ; j < M ; j += 2) /* step through in pairs and phase if necessary */
	if (x[j] + x[j+1] == 1)	 /* need to phase - add new left score to stored right score */
	  { double score = matchPhaseScore (&match[j], uRef, &sp, k) ;
	    score += arr(rScoreStack[j], --arrayMax(rScoreStack[j]), double) ;
	    if (score > 0) { x[j] = 1 ; x[j+1] = 0 ; } else { x[j] = 0 ; x[j+1] = 1 ; }
	  }
      matchUpdate (match, x, M, uRef, k) ;
      pbwtCursorForwardsRead (uOld) ;
      for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ; /* write into pNew */
      pbwtCursorWriteForwards (uNew) ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  pbwtCursorDestroy (uRef) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  for (j = 0 ; j < M ; j += 2)
    { arrayDestroy (qScoreStack[j]) ; arrayDestroy (rScoreStack[j]) ; }
  free (match) ; free (x) ; free (qScoreStack) ; free (rScoreStack) ;
  return pNew ;
}

/****************** phase using maximal matches at homozygous sites  ****************/

typedef struct { int jRef ; int start ; int end ; } MatchSegment ;
/* matches are semi-open [start,end) so length is end-start */

static Array *maxMatch = 0 ;	/* of MatchSegment */

static void reportMatch (int iq, int jRef, int start, int end)
{
  MatchSegment *ms = arrayp (maxMatch[iq], arrayMax(maxMatch[iq]), MatchSegment) ;
  ms->jRef = jRef ; ms->start = start ; ms->end = end ;
}

static PBWT *referencePhase2 (PBWT *pOld, PBWT *pRef)
{
  /* Strategy of this approach is to make two passes through the data.
     The first one will make for each input sequence in pOld a set of maximal 
     match segments in pRef at sites that are homozygous in the input sequence.
     The second evaluates whether to flip consecutive sites or not at het sites
     guided by a weighted sum of information from overlapping match segments.
     This method is O(pRef->N * pRef->M * pOld->M) time and O(pRef->M * pOld->M) space,
     rather than previous methods that are O(pRef->N * (pRef->M + pOld->M)) time
     I think, though with a worse constant.
  */
  int i, j, k ;
  uchar *xOld = myalloc (pOld->M, uchar) ;   /* pOld values in original sort order */
  uchar *xRef = myalloc (pRef->M, uchar) ;   /* pRef values in original sort order */
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  maxMatch = myalloc (pOld->M, Array) ; /* remember maxMatch is a package global */
  for (j = 0 ; j < pOld->M ; j+=2) maxMatch[j] = arrayCreate (1024, MatchSegment) ;

  /* first pass: build maxMatch from homologous sites */
  /* use the strategy in matchSequencesSweep() */
  PbwtCursor **uHom = myalloc (pOld->M, PbwtCursor*) ;
  for (j = 0 ; j < pOld->M ; j+=2) uHom[j] = pbwtNakedCursorCreate (pRef->M, 0) ;
  int *f = mycalloc (uOld->M, int) ; /* first location in uHom[j] of longest match to j'th query */
  int *d = mycalloc (uOld->M, int) ; /* start of longest match to j'th query */
  for (k = 0 ; k < pOld->N ; ++k)
    { for (j = 0 ; j < pOld->M ; ++j) xOld[uOld->a[j]] = uOld->y[j] ;
      for (j = 0 ; j < pRef->M ; ++j) xRef[uRef->a[j]] = uRef->y[j] ;
      for (j = 0 ; j < pOld->M ; j += 2)
	{ int xx ;
	  if (xOld[j] == xOld[j+1]) /* homozygote */
	    { int xx = xOld[j] ;
	      for (i = 0 ; i < pRef->M ; ++i) uHom[j]->y[i] = xRef[uHom[j]->a[i]] ;
	      if (uHom[j]->y[f[j]] != xx) /* ie match does not extend */
		{ int iPlus = f[j] ;
		  while (++iPlus < pRef->M && uHom[j]->d[iPlus] <= d[j])
		    if (uHom[j]->y[iPlus] == xx) { f[j] = iPlus ; goto DONE ; }
		  /* if not, then report these matches */
		  for (i = f[j] ; i < iPlus ; ++i) reportMatch (j, uHom[j]->a[i], d[j], k) ;
		  /* then find new top longest match that can be extended */
		  /* extend out interval [iMinus, iPlus] until we find this best match */
		  int iMinus = f[j] ; /* an index into *uHom[j] less than f[j] */
		  int dPlus = (iPlus < pRef->M) ? uHom[j]->d[iPlus] : k+1 ;
		  int dMinus = uHom[j]->d[iMinus] ;
		  while (TRUE)
		    if (dMinus <= dPlus)
		      { i = -1 ;	/* impossible value */
			while (uHom[j]->d[iMinus] <= dMinus) /* d[0]=k+1 stops underflow */
			  if (uHom[j]->y[--iMinus] == xx) i = iMinus ;
			if (i >= 0) { f[j] = i ; d[j] = dMinus ; goto DONE ; }
			dMinus = uHom[j]->d[iMinus] ;
		      }
		    else		/* dPlus < dMinus */
		      { while (iPlus < pRef->M && uHom[j]->d[iPlus] <= dPlus)
			  if (uHom[j]->y[iPlus] == xx) { f[j] = iPlus ; d[j] = dPlus ; goto DONE ; }
			  else ++iPlus ;
			dPlus = (iPlus == uRef->M) ? k : uHom[j]->d[iPlus] ;
			if (!iMinus && iPlus == uRef->M) 
			  die ("no match to query %d value %d at site %d", j, xx, k) ;
		      }
		}
	    DONE: ;
	      /* next move forwards cursor and update match location */
	      int cc = uHom[j]->c ; pbwtCursorForwardsADU (uHom[j], k) ;
	      f[j] = xx ? cc + (f[j] - uHom[j]->u[f[j]]) : uHom[j]->u[f[j]] ;
	    }
	}
      pbwtCursorForwardsRead (uOld) ;
      pbwtCursorForwardsRead (uRef) ;
    }

  pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uRef) ;
  for (j = 0 ; j < pOld->M ; j+=2) free (uHom[j]) ; free (uHom) ;

  fprintf (stderr, "phaseReference homozygote match pass complete: ") ; timeUpdate() ;

  /* now second pass: phase using matches at homozygous sites
     and encode the rephased haps back into pNew */
  PBWT *pNew = pbwtCreate (pOld->M, pOld->N) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  uchar **xRefStore = mycalloc (pRef->N, uchar*) ; /* store of xRef columns if Count > 0 */
  int *xRefCount = mycalloc (pRef->N, int) ;       /* count of active links to Store */
  int *kLast = myalloc (pOld->M, int) ;     /* k position of last het */
  for (j = 0 ; j < pOld->M ; ++j) kLast[j] = -1 ;
  uchar *xLast = mycalloc (pOld->M, uchar) ; /* for even j, the value at the last het */
  int *firstSeg = mycalloc (pOld->M, int) ;
  for (k = 0 ; k < pOld->N ; ++k)
    { for (j = 0 ; j < pOld->M ; ++j) xOld[uOld->a[j]] = uOld->y[j] ;
      for (j = 0 ; j < pRef->M ; ++j) xRef[uRef->a[j]] = uRef->y[j] ;
      for (j = 0 ; j < pOld->M ; j += 2)
	if (xOld[j] != xOld[j+1])
	  { if (kLast[j] > -1)	/* phase it unless it is the first het */
	      { double bit, sum = 0, score = 0 ;
		uchar *xRefLast = xRefStore[kLast[j]] ;
		/* loop over match segments, making weighted sum of how many flip */
		MatchSegment *m = arrp(maxMatch[j],firstSeg[j],MatchSegment) ;
		MatchSegment *mStop = arrp(maxMatch[j],arrayMax(maxMatch[j]),MatchSegment) ;
		while (m->end <= k && m < mStop) { ++m ; ++firstSeg[j] ; }
		while (m->start <= k && m < mStop)
		  { bit = (k - m->start + 1) * (m->end - k) ; if (!bit) bit = 1 ;
		    sum += bit ;
		    if (xRef[m->jRef] != xRefLast[m->jRef]) score += bit ;
		    ++m ;
		  }
		/* here is the actual phasing */
		if (score/sum > 0.5)
		  { xOld[j+1] = xLast[j] ; xOld[j] = 1 - xOld[j+1] ; xLast[j] = xOld[j] ; }
		else
		  { xOld[j] = xLast[j] ; xOld[j+1] = 1 - xOld[j] ; }
		/* decrement the stored Ref counter and free if necessary */
		if (!--xRefCount[kLast[j]]) free (xRefStore[kLast[j]]) ;
	      }
	    kLast[j] = k ;
	    /* increment the stored Ref counter, and store if necessary */
	    if (!xRefCount[k]++) 
	      { xRefStore[k] = myalloc (pRef->M, uchar) ; 
		memcpy (xRefStore[k], xRef, pRef->M) ;
	      }
	  }
      for (j = 0 ; j < pOld->M ; ++j) uNew->y[j] = xOld[uNew->a[j]] ;
      pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uOld) ;
      pbwtCursorForwardsRead (uRef) ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  free (xOld) ; free (xRef) ;
  pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uRef) ; pbwtCursorDestroy (uNew) ;
  for (j = 0 ; j < pOld->M ; j+=2) free (maxMatch[j]) ; free (maxMatch) ;
  for (k = 0 ; k < pRef->N ; ++k) if (xRefCount[k]) free (xRefStore[k]) ; free (xRefStore) ;
  free (xRefCount) ; free (xLast) ; free (firstSeg) ;
  return pNew ;
}

/************* top level selector for referencePhase *****************/

PBWT *referencePhase (PBWT *pOld, char *fileNameRoot)
{
  fprintf (stderr, "phase against reference %s\n", fileNameRoot) ;
  if (!pOld || !pOld->yz || !pOld->sites) 
    die ("referencePhase called without existing pbwt with sites") ;
  PBWT *pRef = pbwtReadAll (fileNameRoot) ;
  if (!pRef->sites) die ("new pbwt %s in referencePhase has no sites", fileNameRoot) ;
  if (strcmp(pOld->chrom,pRef->chrom))
    die ("mismatching chrom in referencePhase: old %s, new %s", pRef->chrom, pOld->chrom) ;

  /* reduce both down to the intersecting sites */
  pOld = pbwtSelectSites (pOld, pRef->sites, FALSE) ;
  pRef = pbwtSelectSites (pRef, pOld->sites, FALSE) ;
  if (!pOld->N) die ("no overlapping sites in referencePhase") ;

  fprintf (stderr, "Phase preliminaries: ") ; timeUpdate() ;

  PBWT *pNew = referencePhase2 (pOld, pRef) ;
  fprintf (stderr, "Phasing complete: ") ; timeUpdate() ;
  fprintf (stderr, "After phasing: ") ; phaseCompare (pNew, pOld) ;

  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->sites = pOld->sites ; pOld->sites = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtDestroy (pRef) ;
  return pNew ;
}

/*********************************************************************************************/
/**** standard genotype imputation - based heavily on referencePhase1 - should rationalise ***/

static void matchUpdateNext (MatchInfo **pMatch, uchar *xx, int M, PbwtCursor *u, int k)
/* same as matchUpdate but for imputation we store the new info in the next array of MatchInfo */  
{
  int i, j ; MatchInfo *mOld, *mNew ; uchar *x ;

  for (j = 0, mOld = *pMatch, mNew = *(pMatch+1), x = xx ; j < M ; ++j, ++mOld, ++mNew, ++x)
    { mNew->dplus = mOld->dplus ;
      for (i = mOld->i ; i < u->M ; ++i)
	{ if (u->y[i] == *x) break ;
	  if (u->d[i] > mNew->dplus) mNew->dplus = u->d[i] ;
	}
      mNew->dminus = mOld->dminus ;
      for (i = mOld->i-1 ; i >= 0 ; --i)
	{ if (u->y[i] == *x) break ;
	  if (u->d[i] > mNew->dminus) mNew->dminus = u->d[i] ;
	}
    }
  int c = u->c ; pbwtCursorForwardsReadADU (u, k) ;
  for (j = 0, mOld = *pMatch, mNew = *(pMatch+1), x = xx ; j < M ; ++j, ++mOld, ++mNew, ++x)
    mNew->i = *x ? c + mOld->i - u->u[mOld->i] : u->u[mOld->i] ;
}

static inline double matchImputeScore (MatchInfo *m, PbwtCursor *uRef, PbwtCursor *uFrame,
				       int *aInv, ScoreParams *sp, int k)
{
  double score = 0 ;
  if (m->i > 0) 
    score += (2*uRef->y[aInv[uFrame->a[m->i-1]]] - 1) * (sp->alpha + sp->beta*(k - m->dminus)) ;
  if (m->i < uRef->M) 
    score += (2*uRef->y[aInv[uFrame->a[m->i]]] - 1) * (sp->alpha + sp->beta * (k - m->dplus)) ;
  return score ;
}

static PBWT *referenceImpute1 (PBWT *pOld, PBWT *pRef, PBWT *pFrame)
{
/* Strategy: first pass back through the pbwts to build imputation information from right.
   Then pass forwards building information from left and imputing from both.
   The problem here is the potential size of the cache of right imputation information.
   Full data size is N (length) * pOld->M.  Could be more sophisticated.
*/

  fprintf (stderr, "REFERENCE IMPUTE 1: ") ;

  /* declarations and initialisations */
  int M = pOld->M, N = pOld->N ;
  int i, j, k ;		    /* indices: ref/frame, new sample, site */
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, FALSE) ;   /* cursor on old - run backwards on forwards pbwt (we just want the values, not the history, and may not have reverse) */
  PbwtCursor *uFrameR = pbwtCursorCreate (pFrame, FALSE, TRUE) ;   /* cursor on frame - run forwards on reverse pbwt*/
  MatchInfo **matchRstore = myalloc (N+1, MatchInfo*) ; /* impute info from right per site */
  for (k = 0 ; k <= N ; ++k) matchRstore[k] = mycalloc (M, MatchInfo) ;
  MatchInfo **matchR = matchRstore ;
  for (j = 0 ; j < M ; ++j)	    /* initialise first matchR[].i randomly */
    (*matchR)[j].i = rand() * (double)pFrame->M / RAND_MAX ;
  uchar *x = myalloc (M, uchar) ;   /* current uOld values in original sort order */

  /* initial reverse pass */
  for (k = 0 ; k < N ; ++k)
    { pbwtCursorReadBackwards (uOld) ;
      for (j = 0 ; j < M ; ++j)	x[uOld->a[j]] = uOld->y[j] ; /* extract genotype */
      /* Update match, but unlike in referencePhase we store the new match in next matchR */
      matchUpdateNext (matchR++, x, M, uFrameR, k) ; /* includes cursorForwards (uFrameR) */
    }
  fprintf (stderr, "reverse pass complete\n") ;
  if (isCheck) { fprintf (stderr, "After reverse pass: ") ; timeUpdate() ; }

 /* Now the forward pass, building pNew */
  PbwtCursor *uFrameL = pbwtCursorCreate (pFrame, TRUE, TRUE) ; /* now run forwards on forwards pbwt */
  MatchInfo *matchL = mycalloc (M, MatchInfo) ;
  for (j = 0 ; j < M ; ++j) matchL[j].i = (*matchR)[j].i ; /* initialise from final matchR */
  PBWT *pNew = pbwtCreate (M, pRef->N) ;	/* this will hold the imputed sequence */
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  PbwtCursor *uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  int *aRefInv = myalloc (M, int) ;  /* holds the inverse mapping from uRef->a[i] -> i */
  double score ; ScoreParams sp ;    /* will need for scoring */

  /* first impute all sites to the left of the first frame site, from the last right match */
  int xSiteNextFrame = arrp(pFrame->sites,0,Site)->x ;
  int kRef = 0 ; 
  while (arrp(pRef->sites,kRef,Site)->x < xSiteNextFrame)
    { sp = scoreParams(uRef) ;
      for (i = 0 ; i < pRef->M ; ++i) aRefInv[uRef->a[i]] = i ;
      for (j = 0 ; j < M ; ++j)
	{ double score = matchImputeScore (&((*matchR)[j]), uRef, uFrameR, aRefInv, &sp, N) ;
	  x[j] = score > 0 ? 1 : 0 ;
	}
      for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ; /* write into pNew */
      pbwtCursorWriteForwards (uNew) ;
      arrp(pRef->sites,kRef,Site)->freq = (uRef->M - uRef->c) / (double) pRef->M ;
      pbwtCursorForwardsRead (uRef) ; ++kRef ;
    }

  int kFrame = 0 ; 
  while (kFrame < N)
    { if (xSiteNextFrame != arrp(pRef->sites, kRef, Site)->x) die ("error in refImpute()") ;

      /* first fill in the value at the frame site */
      for (j = 0 ; j < M ; ++j) x[uOld->a[j]] = uOld->y[j] ; /* extract genotype */
      for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ; pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uRef) ; ++kRef ;

      /* move forwards the pbwt and match info across the frame site */
      pbwtCursorForwardsRead (uOld) ;
      matchUpdate (matchL, x, M, uFrameL, kFrame) ;	 /* update left match infos */
      pbwtCursorForwardsRead (uFrameL) ; ++kFrame ;
      --matchR ;					 /* go back to previous right match */
      pbwtCursorReadBackwards (uFrameR) ;

      if (kFrame == N) break ;

      /* now impute sites between this frame site and the next */
      xSiteNextFrame = arrp(pFrame->sites,kFrame,Site)->x ;
      while (arrp(pRef->sites,kRef,Site)->x < xSiteNextFrame)
	{ sp = scoreParams(uRef) ;
	  for (i = 0 ; i < pRef->M ; ++i) aRefInv[uRef->a[i]] = i ;
	  for (j = 0 ; j < M ; ++j)
	    { double score = matchImputeScore (matchL, uRef, uFrameL, aRefInv, &sp, kFrame-1) ;
	      score += matchImputeScore (&((*matchR)[j]), uRef, uFrameR, aRefInv, &sp, N-kFrame) ;
	      x[j] = score > 0 ? 1 : 0 ;
	    }
	  for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ; /* write into pNew */
	  pbwtCursorWriteForwards (uNew) ;
	  arrp(pRef->sites,kRef,Site)->freq =  (uRef->M - uRef->c) / (double) pRef->M ;
	  pbwtCursorForwardsRead (uRef) ; ++kRef ;
	}
    }

  /* finally impute from the last frame site to the end */
  while (kRef < pRef->N)
    { sp = scoreParams(uRef) ;
      for (i = 0 ; i < pRef->M ; ++i) aRefInv[uRef->a[i]] = i ;
      for (j = 0 ; j < M ; ++j)
	{ double score = matchImputeScore (matchL, uRef, uFrameL, aRefInv, &sp, N-1) ;
	  x[j] = score > 0 ? 1 : 0 ;
	}
      for (j = 0 ; j < M ; ++j) uNew->y[j] = x[uNew->a[j]] ; /* write into pNew */
      pbwtCursorWriteForwards (uNew) ;
      arrp(pRef->sites,kRef,Site)->freq =  (uRef->M - uRef->c) / (double) pRef->M ;
      pbwtCursorForwardsRead (uRef) ; ++kRef ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  pbwtCursorDestroy (uFrameL) ; pbwtCursorDestroy (uFrameR) ;
  pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ; pbwtCursorDestroy (uRef) ;
  free (aRefInv) ; free (matchL) ; 
  for (k = 0 ; k < N ; ++k) free (matchRstore[k]) ; free (matchRstore) ;
  return pNew ;
}

/********** impute using maximal matches - use matchSequencesSweep() to find them *********/

static PBWT *referenceImpute2 (PBWT *pOld, PBWT *pRef, PBWT *pFrame)
/* require pOld and pFrame to have the same sites, and these are a subset of sites of pRef */
{
  int i, j, k ;

  fprintf (stderr, "Reference impute using maximal matches: ") ;

  /* build the array of maximal matches into pFrame for each sequence in pOld */
  maxMatch = myalloc (pOld->M, Array) ;
  for (j = 0 ; j < pOld->M ; ++j) maxMatch[j] = arrayCreate (1024, MatchSegment) ;
  matchSequencesSweep (pFrame, pOld, reportMatch) ;
  for (j = 0 ; j < pOld->M ; ++j)		    /* add terminating element to arrays */
    { MatchSegment *ms = arrayp(maxMatch[j],arrayMax(maxMatch[j]),MatchSegment) ;
      ms->jRef = ms[-1].jRef ; ms->end = pOld->N+1 ; ms->start = pOld->N ;
    }

  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  PBWT *pNew = pbwtCreate (pOld->M, pRef->N) ;	/* this will hold the imputed sequence */
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  uchar *x = myalloc (pOld->M, uchar) ;   /* current uOld values in original sort order */
  int *aRefInv = myalloc (pRef->M, int) ;  /* holds the inverse mapping from uRef->a[i] -> i */
  int *firstSeg = mycalloc (pOld->M, int) ; /* position in maxMatch to start looking at */

  int kOld = 0, kRef = 0 ;
  while (kRef < pRef->N)
    { if (arrp(pRef->sites,kRef,Site)->x == arrp(pFrame->sites,kOld,Site)->x)
	{ pbwtCursorForwardsRead (uOld) ; ++kOld ;
	  for (j = 0 ; j < pOld->M ; ++j)
	    while (kOld >= arrp(maxMatch[j],firstSeg[j],MatchSegment)->end) ++firstSeg[j] ;
	}
      for (i = 0 ; i < pRef->M ; ++i) aRefInv[uRef->a[i]] = i ;
      for (j = 0 ; j < pOld->M ; ++j)
	{ /* need to look at overlapping matches */
	  double bit = 0, sum = 0, score = 0 ;
	  MatchSegment *m = arrp(maxMatch[j],firstSeg[j],MatchSegment) ;
	  MatchSegment *mStop = arrp(maxMatch[j],arrayMax(maxMatch[j]),MatchSegment) ;
	  while (m->start <= kOld && m < mStop)
	    { bit = (kOld - m->start + 1) * (m->end - kOld) ; if (!bit) bit = 1 ;
	      if (bit < 0)
		die ("refImpute2 kRef %d kOld %d j %d m->start %d m->end %d firstSeg %d\n", 
		     kRef, kOld, j, m->start, m->end, firstSeg[j]) ;
	      sum += bit ;
	      if (uRef->y[aRefInv[m->jRef]]) score += bit ;
	      ++m ;
	    }
	  if (sum == 0) 
	    { 
	      x[j] = 0 ;
	    }
	  else 
	    x[j] = (score/sum > 0.5) ? 1 : 0 ;
	}
      for (j = 0 ; j < pOld->M ; ++j) uNew->y[j] = x[uNew->a[j]] ; /* transfer to uNew */
      pbwtCursorWriteForwards (uNew) ;
      arrp(pRef->sites,kRef,Site)->freq = (uRef->M - uRef->c) / (double) pRef->M ;
      pbwtCursorForwardsRead (uRef) ; ++kRef ;
      if (isCheck && !(kRef % 10000)) fprintf (stderr, " %d %d\n", kRef, kOld) ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uRef) ; pbwtCursorDestroy (uNew) ;
  free (aRefInv) ; free (firstSeg) ;
  for (j = 0 ; j < pOld->M ; ++j) free (maxMatch[j]) ; free (maxMatch) ;
  return pNew ;
}

/*********************************************************************/

PBWT *referenceImpute (PBWT *pOld, char *fileNameRoot)
{
  /* Preliminaries */
  fprintf (stderr, "impute against reference %s\n", fileNameRoot) ;
  if (!pOld || !pOld->yz || !pOld->sites) 
    die ("referenceImpute called without existing pbwt with sites") ;
  PBWT *pRef = pbwtReadAll (fileNameRoot) ;
  if (!pRef->sites) die ("new pbwt %s in referencePhase has no sites", fileNameRoot) ;
  if (strcmp(pOld->chrom,pRef->chrom))
    die ("mismatching chrom in referenceImpute: old %s, new %s", pRef->chrom, pOld->chrom) ;

  /* identify the intersecting sites */
  pOld = pbwtSelectSites (pOld, pRef->sites, FALSE) ;
  if (!pOld->N) die ("no overlapping sites in referenceImpute") ;
  PBWT *pFrame = pbwtSelectSites (pRef, pOld->sites, TRUE) ; /* keep the full ref to impute to */
  if (pFrame->N == pRef->N)
    { fprintf (stderr, "No additional sites to impute in referenceImpute\n") ;
      pbwtDestroy (pFrame) ; pbwtDestroy (pRef) ;
      return pOld ;
    }
  pbwtBuildReverse (pFrame) ;	/* we need the reverse reference pbwt below */
  if (!pOld->aFend) die ("pOld has no aFend in referenceImpute - your pbwt was made by a previous version of the code; buildReverse and resave the forwards pbwt") ;

  fprintf (stderr, "Imputation preliminaries: ") ; timeUpdate() ;

  PBWT *pNew = referenceImpute2 (pOld, pRef, pFrame) ;
  pNew->sites = pRef->sites ; pRef->sites = 0 ; 
  pNew->chrom = pRef->chrom ; pRef->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtDestroy (pFrame) ; pbwtDestroy (pRef) ;
  return pNew ;
}

/*********************************************************************/

PBWT *imputeMissing (PBWT *p)
{
  return p ;
}

/*********************************************************************/

void genotypeCompare (PBWT *p, char *fileNameRoot)
{
  fprintf (stderr, "compare genotypes to reference %s\n", fileNameRoot) ;
  if (!p || !p->yz || !p->sites) 
    die ("genotypeCompare called without existing pbwt with sites") ;
  PBWT *pRef = pbwtReadAll (fileNameRoot) ;
  if (!pRef->sites) die ("new pbwt %s in genotypeCompare has no sites", fileNameRoot) ;
  if (p->M != pRef->M) die ("mismatch of old M %d to ref M %d", p->M, pRef->M) ;
  if (p->N != pRef->N) warn ("mismatch of old N %d to ref N %d", p->N, pRef->N) ;
  
  /* reduce both down to the intersecting sites */
  PBWT *pFrame = pbwtSelectSites (p, pRef->sites, TRUE) ;
  pRef = pbwtSelectSites (pRef, p->sites, FALSE) ;
  if (!pFrame->N) die ("no overlapping sites in genotypeCompare") ;
  if (strcmp(pFrame->chrom,pRef->chrom)) die ("mismatch chrom %s to ref %f", pFrame->chrom, pRef->chrom) ;

  genotypeComparePbwt (pFrame, pRef) ;
  
  pbwtDestroy (pRef) ; pbwtDestroy(pFrame) ;
}

static void genotypeComparePbwt (PBWT *p, PBWT *q)
{
  int i, j, k, ff ;
  long n[16][9] ; for (i = 16 ; i-- ;) for (j = 9 ; j-- ;) n[i][j] = 0 ;
  long *ns = mycalloc (9*p->M, long) ;

  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  PbwtCursor *uq = pbwtCursorCreate (q, TRUE, TRUE) ;
  uchar *xp = myalloc (p->M, uchar), *xq = myalloc (q->M, uchar) ;
  for (k = 0 ; k < p->N ; ++k)
    { double f = (p->M - up->c) / (double)p->M ; 
      if (arrp (p->sites,k,Site)->freq) f = arrp (p->sites,k,Site)->freq ;
      for (ff = 0 ; f*100 > fBound[ff] ; ++ff) ;
      for (j = 0 ; j < p->M ; ++j) { xp[up->a[j]] = up->y[j] ; xq[uq->a[j]] = uq->y[j] ; }
      for (j = 0 ; j < p->M ; j += 2) 
	{ i = 3*(xp[j]+xp[j+1]) + xq[j]+xq[j+1] ;
	  ++n[ff][i] ;
	  ++ns[9*j + i] ;
	}
      pbwtCursorForwardsRead (up) ; pbwtCursorForwardsRead (uq) ;
    }
  pbwtCursorDestroy (up) ; pbwtCursorDestroy (uq) ;

  /* report */
  for (ff = 0 ; ff < 16 ; ++ff)
    { printf ("%-5.1f", fBound[ff]) ;
      double tot = 0, xbar, ybar, r2, x2, y2 ;
      for (i = 0 ; i < 9 ; ++i) { printf ("\t%ld ", n[ff][i]) ; tot += n[ff][i] ; }
      if (tot) 
	{ xbar = (n[ff][3] + n[ff][4] + n[ff][5] + 2*(n[ff][6] + n[ff][7] + n[ff][8])) / tot ;
	  x2 = (n[ff][3] + n[ff][4] + n[ff][5] + 4*(n[ff][6] + n[ff][7] + n[ff][8])) / tot ;
	  ybar = (n[ff][1] + n[ff][4] + n[ff][7] + 2*(n[ff][2] + n[ff][5] + n[ff][8])) / tot ;
	  y2 = (n[ff][1] + n[ff][4] + n[ff][7] + 4*(n[ff][2] + n[ff][5] + n[ff][8])) / tot ;
	  r2 = (n[ff][4] + 2*(n[ff][5] + n[ff][7]) + 4*n[ff][8]) / tot ;
	  r2 = (r2 - xbar*ybar)/sqrt((x2 - xbar*xbar)*(y2 - ybar*ybar)) ;
	  printf ("\tx,y,r2\t%.4f\t%.4f\t%.4f\n", xbar, ybar, r2) ;
	}
      else
	putchar ('\n') ;
    }
  int hist[101] ; for (i = 101 ; i-- ;) hist[i] = 0 ;
  for (j = 0 ; j < p->M ; ++j)
    { double tot = 0, xbar, ybar, r2, x2, y2 ;
      for (i = 0 ; i < 9 ; ++i) tot += ns[9*j+i] ;
      if (tot) 
	{ xbar = (ns[9*j+3] + ns[9*j+4] + ns[9*j+5] + 2*(ns[9*j+6] + ns[9*j+7] + ns[9*j+8])) / tot ;
	  x2 = (ns[9*j+3] + ns[9*j+4] + ns[9*j+5] + 4*(ns[9*j+6] + ns[9*j+7] + ns[9*j+8])) / tot ;
	  ybar = (ns[9*j+1] + ns[9*j+4] + ns[9*j+7] + 2*(ns[9*j+2] + ns[9*j+5] + ns[9*j+8])) / tot ;
	  y2 = (ns[9*j+1] + ns[9*j+4] + ns[9*j+7] + 4*(ns[9*j+2] + ns[9*j+5] + ns[9*j+8])) / tot ;
	  r2 = (ns[9*j+4] + 2*(ns[9*j+5] + ns[9*j+7]) + 4*ns[9*j+8]) / tot ;
	  r2 = (r2 - xbar*ybar)/sqrt((x2 - xbar*xbar)*(y2 - ybar*ybar)) ;
	  if (r2 < 0) r2 = 0 ;
	  ++hist[(int)(100*r2)] ;
	}
    }
  if (hist[100]) printf ("%d samples with r2 == 1.0\n", hist[100]) ;
  for (i = 100 ; i-- ; ) 
    if (hist[i])
      printf ("%d samples with %.2f <= r2 < %.2f\n", hist[i], (i-1)*0.01, i*0.01) ;
}

/*********************************************************************/

void imputeEvaluate (PBWT *p, char *fileNameRoot)
{
}

/*********** routines to corrupt data to explore robustness *************/

PBWT *pbwtCorruptSites (PBWT *pOld, double pSite, double pChange)
{
  int M = pOld->M, N = pOld->N ;
  PBWT *pNew = pbwtCreate (M, N) ;
  int rSite = pSite*RAND_MAX, rChange = pChange*RAND_MAX ;
  double rFac = RAND_MAX / (double) M ;
  int k, i ;
  uchar *x ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  int nChange = 0 ;

  if (!pOld || !pOld->yz) die ("corruptSites without an existing pbwt") ;
  if (pSite <= 0 || pSite > 1 || pChange <= 0 || pChange > 1)
    die ("pSite %f, pChange %f for corruptSites out of range\n", pSite, pChange) ;

  x = myalloc (M, uchar) ;

  for (k = 0 ; k < N ; ++k)
    { for (i = 0 ; i < M ; ++i) x[uOld->a[i]] = uOld->y[i] ;
	  
      for (i = 0 ; i < M ; ++i) uNew->y[i] = x[uNew->a[i]] ;
      if (rand() < rSite)
	for (i = 0 ; i < M ; ++i)
	  if (rand() < rChange)
	    { uchar old = uNew->y[i] ;
	      uNew->y[i] = (rand() < uOld->c*rFac) ? 0 : 1 ;
	      if (old != uNew->y[i])
		++nChange ;
	    }
      pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uOld) ;
    }  
  pbwtCursorToAFend (uNew, pNew) ;

  fprintf (stderr, "corruptSites with pSite %f, pChange %f changes %.4f of values\n", 
	   pSite, pChange, nChange/(N*(double)M)) ;

  pNew->sites = pOld->sites ; pOld->sites = 0 ; 
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ;
  return pNew ;
}

PBWT *pbwtCorruptSamples (PBWT *pOld, double pSample, double pChange)
{
  int M = pOld->M, N = pOld->N ;
  PBWT *pNew = pbwtCreate (M, N) ;
  int rSample = pSample*RAND_MAX, rChange = pChange*RAND_MAX ;
  double rFac = RAND_MAX / (double) M ;
  int k, i ;
  uchar *x ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  BOOL *isCorrupt = myalloc (M, BOOL) ;
  int nChange = 0 ;

  if (!pOld || !pOld->yz) die ("corruptSites without an existing pbwt") ;
  if (pSample <= 0 || pSample > 1 || pChange <= 0 || pChange > 1)
    die ("pSample %f, pChange %f for corruptSites out of range\n", pSample, pChange) ;

  x = myalloc (M, uchar) ;

  for (i = 0 ; i < M ; ++i)
    isCorrupt[i] = (rand() < rSample) ;

  for (k = 0 ; k < N ; ++k)
    { for (i = 0 ; i < M ; ++i) x[uOld->a[i]] = uOld->y[i] ;
	  
      for (i = 0 ; i < M ; ++i) 
	if (isCorrupt[i] && rand() < rChange)
	  { uNew->y[i] = (rand() < uOld->c*rFac) ? 0 : 1 ;
	    if (uNew->y[i] != x[uNew->a[i]]) ++nChange ;
	  }
	else
	  uNew->y[i] = x[uNew->a[i]] ;
      pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uOld) ;
    }  
  pbwtCursorToAFend (uNew, pNew) ;

  fprintf (stderr, "corruptSamples with pSample %f, pChange %f changes %.4f of values\n",
	   pSample, pChange, nChange/(N*(double)M)) ;
  
  pNew->sites = pOld->sites ; pOld->sites = 0 ; 
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ; free (isCorrupt) ;
  return pNew ;
}

/******************* end of file *******************/

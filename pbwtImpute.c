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
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 22 18:41 2013 (rd)
 * Created: Thu Apr  4 12:02:56 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

/************************* IMPUTATION AND PHASING *************************/

#include <math.h>

void imputeExplore (PBWT *p, int test)
{
  int i, j, k, M = p->M, N = p->N, ff ;
  uchar *y ;
  double xbar, ybar, tot, r2 ;
  static double fBound[] = {0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20, 50, 100.01} ;
  double f ;
  typedef struct {
    long n00, n01 ;		/* neither neighbour is 1, truth is 0 or 1 */
    long n10a, n10b, n11a, n11b ;/* one neighbour is 1; a if d is lower for the 0 neighbour */
    long n20, n21 ;		/* both neighbours are 1 */
  } TestStat ;
  TestStat *testStat = mycalloc (15, TestStat), *t ;
  typedef long Counts[4] ;
  Array dHist = arrayCreate (1000, Counts) ;
  Counts cSimple, cCond0, cCond1 ;
  int *n0 = myalloc (M, int), *n1 = myalloc (M, int) ;
  uchar *x = myalloc (M, uchar) ;
  static long c0[15][5], c1[15][5] ;

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
    for (j = 0 ; j < 15 ; ++j)
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
    for (j = 0 ; j < 15 ; ++j)
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

  q = pbwtCreate (M) ; q->yz = arrayCreate (arrayMax(p->yz), uchar) ; q->N = N ;

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
	{ PBWT *r = pbwtCreate (M) ; r->N = N ;
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

PBWT *referencePhase (PBWT *pOld, char *fileNameRoot)
{
  if (!pOld || !pOld->yz || !pOld->sites) 
    die ("referencePhase called without existing pbwt with sites") ;
  PBWT *pNew = pbwtReadAll (fileNameRoot) ;
  if (!pNew->sites) die ("new pbwt %s in referencePhase has no sites", fileNameRoot) ;
  if (strcmp(pNew->chrom,pOld->chrom))
    die ("mismatching chom in referencePhase: old %s, new %s", pOld->chrom, pNew->chrom) ;

  /* reduce both down to the intersecting sites */
  pOld = pbwtSelectSites (pOld, pNew->sites) ;
  if (!pOld->zz) pbwtBuildReverse (pOld) ; /* need reverse below */
  pNew = pbwtSelectSites (pNew, pOld->sites) ;

  /* now phase the new sites against the old sites */
  PbwtCursor *ufOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;   /* forwards cursor on old */
  PbwtCursor *ufNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;   /* forwards cursor on new */
  PBWT *q = pbwtCreate (pNew->M) ; q->N = pNew->N ;	      /* will hold the new phasing */
  q->yz = arrayCreate (arrayMax(pNew->yz), uchar) ;
  PbwtCursor *uq = pbwtCursorCreate (q, TRUE, TRUE) ;
  int *iq = myalloc (pNew->M, int) ; /* the position in y i'th new hap follow */
  int *dq = mycalloc (pNew->M, int) ; /* start of match to following position */
  uchar *x = myalloc (pNew->M, uchar) ;
  uchar *yq = myalloc (pNew->M, uchar) ;

  int i, j, k ;			/* site, new sample, old sample */
  for (j = 0 ; j < pNew->M ; ++j) iq[j] = rand() * (double)pNew->M / RAND_MAX ;
  for (k = 0 ; k < pOld->N ; ++k) /* loop forwards */
    { for (j = 0 ; j < pNew->M ; ++j) x[ufNew->a[j]] = ufNew->y[j] ; /* extract genotype */
      for (j = 0 ; j < pNew->M ; j += 2) /* step through in pairs and phase if necessary */
	if (x[j] + x[j+1] == 1)	 /* need to phase */
	  { 
	  }
      int c = ufOld->c ;
      for (j = 0 ; j < pNew->M ; ++j)
	if (iq[j] < pNew->M) yq[j] = ufOld->y[iq[j]] ; else yq[j] = 2 ;
      pbwtCursorForwardsReadADU (ufOld, k) ; /* need the new ufOld->u to update iq[]  */
      for (j = 0 ; j < pNew->M ; ++j)
	{ iq[j] = x[j] ? c + iq[j] - ufOld->u[iq[j]] : ufOld->u[iq[j]] ;
/* Let's try to have iq between 0 and M inclusive be the position just AFTER the query.
   ufOld->u[k] is the number of 0s in y[] before position k.  
   So if x == 0 and there are no 0s before k = iq then iq is mapped to 0 - good.
   If x == 1 and there are no 1s then u[iq] = iq so iq is mapped to c - good.
   If x == 1 and all 1s are before iq then iq is mapped to M - good.
   Mapping of M is correct, to c if 0 else M if 1, so long as u[M] is defined as c - yes.
   What about dq[]?
   If x[j] == y[ainv[iq[j]]] then dq[] does not change and we are done.
   Otherwise, if x == 0 and y == 1 then 
 */
	}
      pbwtCursorForwardsRead (ufNew) ;
      pbwtCursorWriteForwards (uq) ;
    }

  q->sites = pNew->sites ; pNew->sites = 0 ;
  q->samples = pNew->samples ; pNew->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtDestroy (pNew) ;
  return q ;
}

/*********** routines to corrupt data to explore robustness *************/

PBWT *pbwtCorruptSites (PBWT *pOld, double pSite, double pChange)
{
  int M = pOld->M, N = pOld->N ;
  PBWT *pNew = pbwtCreate (M) ;
  int rSite = pSite*RAND_MAX, rChange = pChange*RAND_MAX ;
  double rFac = RAND_MAX / (double) M ;
  int k, i ;
  uchar *x ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  pNew->yz = arrayCreate (arrayMax(pOld->yz), uchar) ; pNew->N = N ;
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

  fprintf (stderr, "corruptSites with pSite %f, pChange %f changes %.4f of values\n", 
	   pSite, pChange, nChange/(N*(double)M)) ;

  if (pOld->sites) { pNew->sites = pOld->sites ; pOld->sites = 0 ; }
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ;
  return pNew ;
}

PBWT *pbwtCorruptSamples (PBWT *pOld, double pSample, double pChange)
{
  int M = pOld->M, N = pOld->N ;
  PBWT *pNew = pbwtCreate (M) ;
  int rSample = pSample*RAND_MAX, rChange = pChange*RAND_MAX ;
  double rFac = RAND_MAX / (double) M ;
  int k, i ;
  uchar *x ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  pNew->yz = arrayCreate (arrayMax(pOld->yz), uchar) ; pNew->N = N ;
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

  fprintf (stderr, "corruptSamples with pSample %f, pChange %f changes %.4f of values\n",
	   pSample, pChange, nChange/(N*(double)M)) ;
  
  if (pOld->sites) { pNew->sites = pOld->sites ; pOld->sites = 0 ; }
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ; free (isCorrupt) ;
  return pNew ;
}

/******************* end of file *******************/

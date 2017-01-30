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
 * Last edited: Dec 14 13:52 2015 (rd)
 * * Sep 22 23:10 2014 (rd): move to 64 bit arrays
 * Created: Thu Apr  4 12:02:56 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

static void genotypeComparePbwt (PBWT *p, PBWT *q) ; /* forward declaration */

/************************* IMPUTATION AND PHASING *************************/

#include <math.h>

static double fBound[] = {0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 90, 100.01} ;

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
    double fsum ;
  } TestStat ;
  TestStat *testStat = mycalloc (17, TestStat), *t ;
  typedef long Counts[4] ;
  Array dHist = arrayCreate (1000, Counts) ;
  Counts cSimple, cCond0, cCond1 ;
  int *n0 = myalloc (M, int), *n1 = myalloc (M, int) ;
  uchar *x = myalloc (M, uchar) ;
  static long c0[17][5], c1[17][5] ;

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
	      fprintf (logFile, "forward-backward mismatch at k %d i %d\n", k, i) ;
	}
      if (k > 0.2*N && k < 0.8*N)     /* ignore ends */
	{ f = (M - u->c) / (double)M ; for (ff = 0 ; f*100 > fBound[ff] ; ++ff) ;
	  t = &testStat[ff] ;
	  t->fsum += f ;
	  memset (n0, 0, M*sizeof(int)) ; memset (n1, 0, M*sizeof(int)) ; 
	  for (i = 1 ; i < M-1 ; ++i)
	    { if (u->y[i-1] && u->y[i+1])
		if (u->y[i]) ++t->n21 ; else ++t->n20 ;
	      else if (!u->y[i-1] && !u->y[i+1])
		if (u->y[i]) ++t->n01 ; else ++t->n00 ;
	      else if ((!u->y[i-1] && u->d[i] < u->d[i+1]) || 
		       (!u->y[i+1] && u->d[i+1] < u->d[i]))
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
	      { if (x[i]) ++c1[ff][n1[i]] ; else ++c0[ff][n1[i]] ; }
	}
      pbwtCursorForwardsReadAD (u, k) ;
    }

  if (test == 1)
    for (j = 0 ; j < 17 ; ++j)
      { t = &testStat[j] ;
	tot = t->n00 + t->n01 + t->n10a + t->n11a + t->n10b + t->n11b + t->n20 + t->n21 ;
	printf ("%-5.1f\t%-7.3f\t00,01\t%ld\t%ld\t10a,11a\t%ld\t%ld\t10b,11b\t%ld\t%ld\t20,21\t%ld\t%ld",
		fBound[j], tot ? ((double)testStat[j].fsum) / tot : 0.0, 
		t->n00, t->n01, t->n10a, t->n11a, t->n10b, t->n11b, t->n20, t->n21) ;
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
    for (j = 0 ; j < 17 ; ++j)
      { printf ("%-5.1f", fBound[j]) ;
	tot = 0 ; 
	for (i = 0 ; i < 5 ; ++i) { tot += c0[j][i] + c1[j][i] ; }
	printf ("\t%-7.3f", tot ? ((double) testStat[j].fsum) / tot : 0.0) ;
	xbar = 0 ; r2 = 0 ;
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
	  if (isCheck && (xp[i]+xp[i+1] != xq[i]+xq[i+1])) 
	    { fprintf (logFile, "phaseCompare mismatch k %d sequence %d\np", k, i) ;
	      uchar **pHaps = pbwtHaplotypes(p), **qHaps = pbwtHaplotypes(q) ;
	      int kk ; 
	      for (kk = 0 ; kk < 40 ; ++kk) 
		fprintf (logFile, " %d", pHaps[i][kk] + pHaps[i+1][kk]) ;
	      fprintf (logFile, "\nq") ;
	      for (kk = 0 ; kk < 40 ; ++kk) 
		fprintf (logFile, " %d", qHaps[i][kk] + qHaps[i+1][kk]) ;
	      fprintf (logFile, "\n") ;
	      die ("phaseCompare mismatch: k %d, i %d, xp %d|%d, xq %d|%d", 
		   k, i, xp[i], xp[i+1], xq[i], xq[i+1]) ;
	    }
	}
      pbwtCursorForwardsRead (up) ;
      pbwtCursorForwardsRead (uq) ;
    }

  fprintf (logFile, "%.1f switches per sample, %.3f per het, %.1f nSwitch1, %.1f nSwitch5\n", 
	   mFac*nSwitch, nSwitch/(double)nHet, mFac*nSwitch1, mFac*nSwitch5) ;

  if (isStats)
    { for (i = 0 ; i < M/2 ; ++i)
	{ printf ("SAMPLE-SWITCH\t%d\t%d", i, nSwitchSample[i]) ;
	  if (p->samples)
	    printf ("\t%s", sampleName(pbwtSample (p, 2*i))) ;
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

/******* phase() is phaseOld() rewritten with new API and method 3 only ******/

PBWT *phaseSweep (PBWT *p, PBWT *ref, BOOL isStart, PBWT *r, int nSparse)
{
  int    i, j, k, kk, kr = 0 ;
  int    M = p->M, N = p->N ;

  if (ref && p->M > ref->M) die ("phaseSweep requires ref->M >= p->M") ;

  /* initialisation */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, isStart) ;
  PBWT *q = pbwtCreate (M, N) ; /* new pbwt */
  PbwtCursor *ur, *uref ;
  if (ref) uref = pbwtCursorCreate (ref, TRUE, isStart) ; 
  if (r) 
    { ur = pbwtCursorCreate (r, TRUE, FALSE) ;
      memcpy (ur->b, r->aRend, M*sizeof(int)) ; /* recover stored locations */
      // for (i = 0 ; i < M ; ++i) ur->b[ur->a[i]] = i ; /* store inverse a in b */
      memcpy (q->aFstart, r->aFend, M*sizeof(int)) ; /* prime uq with final ur */
    }
  PbwtCursor *uq = pbwtCursorCreate (q, TRUE, TRUE) ; 
  for (i = 0 ; i < M ; ++i) uq->b[uq->a[i]] = i ; /* store inverse a in b */
  PbwtCursor **uqq = myalloc (nSparse, PbwtCursor*) ;
  for (kk = 0 ; kk < nSparse ; ++kk)
    { uqq[kk] = pbwtNakedCursorCreate (M, 0) ; 
      for (i = 0 ; i < M ; ++i) uqq[kk]->b[uqq[kk]->a[i]] = i ;
    }

  /* now loop through p phasing into q */
  uchar  *x = myalloc (M, uchar) ;  /* actual haplotypes in original order, from p */
  double *xp = myalloc(M, double) ; /* 2*p(x=1)-1, so 1 if x=1, -1 if x=0, 0 if unknown */
  for (k = 0 ; k < N ; k++)
    { if (!isStart) pbwtCursorReadBackwards (up) ;
      for (i = 0 ; i < M ; ++i) x[up->a[i]] = up->y[i] ;  /* build x from up->y */
      if (isStart) pbwtCursorForwardsRead (up) ; 
      for (i = 0 ; i < M ; ++i) xp[i] = x[i] ? 1.0 : -1.0 ;
      int n2 = 0 ;
      for (i = 0 ; i < M ; i += 2) /* go through x in pairs */
	if (x[i] != x[i+1]) { ++n2 ; xp[i] = xp[i+1] = 0.0 ; }  /* a het */
      double s, thresh = ref ? 0.5 : 2*(nSparse + (r?2:1)) + 0.5 ;
      while (n2 && thresh > 1.0)
	{ int n2Old = n2 ; n2 = 0 ;
	  for (i = 0 ; i < M ; i += 2) /* loop over genotype pairs in original order */
	    if (!xp[i]) /* a het to phase */
	      { s = score0 (uq, xp, i) - score0 (uq, xp, i+1) ;
		if (r) s += score0 (ur, xp,i) - score0 (ur, xp, i+1) ;
		for (kk = 0 ; kk < nSparse ; ++kk) 
		  s += score0 (uqq[kk], xp, i) - score0 (uqq[kk], xp, i+1) ;
		if (s > thresh)  { xp[i] = 1 ; xp[i+1] = -1 ; }
		else if (s < -thresh) { xp[i] = -1 ; xp[i+1] = 1 ; }
		else ++n2 ;
	      }
	  if (n2 == n2Old) { thresh -= 1.0 ; n2Old = M+1 ; }
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

      for (i = 0 ; i < M ; ++i)	x[i] = (xp[i] > 0.0) ? 1 : 0 ;
      for (i = 0 ; i < M ; ++i)	uq->y[i] = x[uq->a[i]] ;
      pbwtCursorWriteForwardsAD (uq, k) ; 
      for (i = 0 ; i < M ; ++i) uq->b[uq->a[i]] = i ;
      kk = k % nSparse ;	/* which of the sparse pbwts to update this time */
      for (i = 0 ; i < M ; ++i) uqq[kk]->y[i] = x[uqq[kk]->a[i]] ;
      pbwtCursorForwardsAD (uqq[kk], k/nSparse) ; 
      for (i = 0 ; i < M ; ++i) uqq[kk]->b[uqq[kk]->a[i]] = i ;
      if (r) 
	{ pbwtCursorReadBackwards (ur) ; 
	  for (i = 0 ; i < M ; ++i) ur->b[ur->a[i]] = i ;
	}
    }
  pbwtCursorToAFend (uq, q) ;
  /* cache uq->b in aRend so we can retrieve from reverse on forwards pass */
  q->aRend = myalloc (q->M, int) ; 
  for (i = 0 ; i < q->M ; ++i) q->aRend[i] = uq->b[i] ;
  /* clean up memory allocated */
  free (x) ; free (xp) ;
  pbwtCursorDestroy (up) ; pbwtCursorDestroy (uq) ; if (r) pbwtCursorDestroy (ur) ;
  for (kk = 0 ; kk < nSparse ; ++kk) pbwtCursorDestroy (uqq[kk]) ; free (uqq) ;
  return q ;
}

PBWT *phase (PBWT *p, int nSparse) /* return rephased p */
{
  if (p->M % 2) die ("phase requires that M = %d is even", p->M) ;
  if (nSparse < 2) nSparse = 2 ;

  phaseInit (p->N) ;		/* sets up some lookup tables */

  if (!p->aFend) pbwtBuildReverse (p) ;	/* needed for old PBWT format data */
  PBWT *r = phaseSweep (p, 0, FALSE, 0, 2) ; /* always reverse sweep wth nSparse 2 */
  if (isCheck)		/* flip p->zz round into p->yz and compare to r */
    { Array yzStore = p->yz ; p->yz = p->zz ;
      int *aFstartStore = p->aFstart ; p->aFstart = p->aRstart ;
      fprintf (logFile, "After reverse pass: ") ; phaseCompare (p, r) ;
      p->yz = yzStore ; p->aFstart = aFstartStore ;
    }
  PBWT *q = phaseSweep (p, 0, TRUE, r, nSparse) ;

  /* compare new phasing to original and report switch rates */
  fprintf (logFile, "After forward pass: ") ; phaseCompare (p, q) ;

  pbwtDestroy (p) ;
  return q ;
}

/******* reference phase using same strategy as phase **********/

PBWT *referencePhase0 (PBWT *p, PBWT *pRef)
{
  int nSparse = 2 ;
  phaseInit (pRef->N) ;
  if (!pRef->aFend) pbwtBuildReverse (pRef) ; /* needed for oldPBWT format data */

  PBWT *r = phaseSweep (p, pRef, FALSE, 0, 2) ; /* always reverse with nSparse 2 */
  if (isCheck)		/* flip p->zz round into p->yz and compare to r */
    { if (!p->zz) pbwtBuildReverse (p) ;
      Array yzStore = p->yz ; p->yz = p->zz ;
      int *aFstartStore = p->aFstart ; p->aFstart = p->aRstart ;
      fprintf (logFile, "After reverse pass: ") ; phaseCompare (p, r) ;
      p->yz = yzStore ; p->aFstart = aFstartStore ;
    }
  PBWT *q = phaseSweep (p, pRef, TRUE, r, nSparse) ;
  return q ;
}

/******* phase a new pbwt against the existing one as a reference *******/

typedef struct { int jRef ; int start ; int end ; } MatchSegment ;
/* matches are semi-open [start,end) so length is end-start */

static Array *maxMatch = 0 ;	/* of MatchSegment */

static void reportMatch (int iq, int jRef, int start, int end)
{
  MatchSegment *ms = arrayp (maxMatch[iq], arrayMax(maxMatch[iq]), MatchSegment) ;
  ms->jRef = jRef ; ms->start = start ; ms->end = end ;
}

/**************** referencePhase4 *************************************************/

  /* This one keeps track for each position j in the reference pbwt of the (normalised) 
     probability of generating the data given that the phasing so far sorts at position j,
     and also of the other haplotype sort position given that phasing.  We keep backtrack
     pointers for each j, and prune them back so that only ones connected to the current 
     values at k are kept live.
     As currently implemented this calculates the full likelihood, but gives the most likely
     path contributing to that.  We could adapt it to full Viterbi, or sampling mode.
  */

  /* within this, try different ways to keep the score and extend */

#define EXTEND4

/******** type declarations ********/

typedef unsigned int TraceBackIndex ;	 /* index into traceBackHeap */

/* #define TRACEBACK_DEBUG */

typedef struct 
{ TraceBackIndex back : 28 ;  /* index into traceBackHeap for cell at previous het site */
  unsigned int value : 1 ;  /* 0 or 1 */
  unsigned int count : 3 ;  /* how many cells at next het site point back to this: 0,1 or 2 */
#ifdef TRACEBACK_DEBUG
  int j, k ;
#endif
} TraceBackCell ;	    /* pack this so it fits into a 4-byte word */

typedef struct {
  int j1 ;		    /* the pair index of j0, which is the index of the cell */
  TraceBackIndex back ;
  float s ;
#ifdef EXTEND0
  float sBest ;
#endif
#ifdef EXTEND1
  int dplus0, dminus0, dmax0, dplus1, dminus1, dmax1 ;
#endif
#ifdef EXTEND2
  int j0min, j0max, j1min, j1max ;
#endif
#if defined(EXTEND3) || defined(EXTEND4)
  int dplus0, dminus0, dplus1, dminus1 ;
#endif
 } PhaseCell ;

 /********** first the package to manage trace back space efficiently **********/

 static Array traceBackHeap = 0 ; /* of TraceBackCell */
 static Array traceBackFreeStack = 0 ; /* of TraceBackIndex */

 #define traceBackCell(tb) arrp(traceBackHeap, (tb), TraceBackCell)

 TraceBackIndex traceBackCreate (TraceBackIndex back, unsigned int value, int j, int k)
 { TraceBackIndex tb ;
   TraceBackCell *tc ;
   if (!isCheck && arrayMax (traceBackFreeStack)) 
     { tb = arr (traceBackFreeStack, --arrayMax(traceBackFreeStack), TraceBackIndex) ;
       tc = arrp (traceBackHeap, tb, TraceBackCell) ;
     }
   else
     { tb = arrayMax(traceBackHeap) ; 
       tc = arrayp (traceBackHeap, tb, TraceBackCell) ; 
     }
   tc->back = back ; tc->value = value ; tc->count = 1 ; 
#ifdef TRACEBACK_DEBUG
   tc->j = j ; tc->k = k ;
#endif
   /* if (isCheck) printf ("traceBackCreate %d at k %d j %d back %d value %d\n", tb, k, j, back, value) ; */
   return tb ;
 }

 inline static void traceBackPrune (TraceBackIndex tb, int k, int j)
 { TraceBackCell *tc = traceBackCell(tb) ;
   while (tb && --tc->count == 0) /* i.e. that was the last active link for that cell */
     { array(traceBackFreeStack, arrayMax(traceBackFreeStack), int) = tb ;
       /*      if (isCheck) printf ("traceBackDestroy %d at k %d j %d\n", tb, k, j) ; */
       tb = tc->back ;
       tc = traceBackCell(tb) ;
     }
 }

 void traceBackInitialise (void)
 { traceBackHeap = arrayReCreate (traceBackHeap, 4096, TraceBackCell) ;
   traceBackFreeStack = arrayReCreate (traceBackFreeStack, 2048, TraceBackIndex) ;
   traceBackCreate (0, 0, 0, 0) ; /* make empty first element, so that later back == 0 only at end*/
 }

#ifdef TRACEBACK_DEBUG
 void traceBackReport (TraceBackIndex tb)
 { while (tb)
     { TraceBackCell *tc = traceBackCell(tb) ;
       printf (" %dk%dj%d", tc->value, tc->k, tc->j) ;
       tb = tc->back ;
     }
 }

 void traceBackCheck (int jq, int k, PhaseCell *pc, int M)	/* beware - very costly */
 {
   int j ;
   int* zCount = mycalloc (arrayMax(traceBackHeap), int) ;
   for (j = 0 ; j < arrayMax(traceBackHeap) ; ++j)
     { TraceBackCell *tc = traceBackCell(j) ;
       if (tc->count) ++zCount[tc->back] ;
     }
   for (j = 0 ; j <= M ; ++j) if (pc[j].s) ++zCount[pc[j].back] ;
   for (j = 1 ; j < arrayMax(traceBackHeap) ; ++j)
     if (zCount[j] && traceBackCell(j)->count != zCount[j])
       { warn ("bad count jq %d, k %d, j %d, zCount[j] %d, traceBackCell(j)->count %d",
	       jq, k, j, zCount[j], traceBackCell(j)->count) ;
	 int jj ;
	 for (jj = 0 ; jj <= M ; ++jj)
	   if (pc[jj].s && pc[jj].back == j) warn (" back[%d] == %d", jj, j) ;
	 for (jj = 0 ; jj < arrayMax(traceBackHeap) ; ++jj)
	   { TraceBackCell *tc = traceBackCell(jj) ;
	     if (tc->count && tc->back == j) 
	       warn (" tc[%d]->back %d, tc[%d]->count %d", jj, j, jj, tc->count) ;

	   }
       }
   free (zCount) ;
 }
#endif	// TRACEBACK_DEBUG

/********** next the calculation of how likely to extend from j with value x **********/

#ifdef EXTEND0
/* first version tries to consider all copying haplotypes with weights proportional to shared length */
/* we pre-calculate contributions from above and below in extendPrepare(), caching in extendScore0/1 */
static char *extendMethodText = "EXTEND0" ;

static long *extendScore0 = 0, *extendScore1 = 0 ;
static int extendScoreSize = 0 ;
static Array extendStack = 0 ;

typedef struct {
  int d ;
  int n0, n1 ;			/* number of y[j] == 0/1 so far with this effective dd */
  int sum0, sum1 ;		/* sum of dd*n up to this element */
} ExtendStruct ;

static ExtendStruct *extendUpdate (int dd, int x)
{
  int i, n0 = 0, n1 = 0 ;	/* n0 and n1 are # 0s,1s with d > dd */
  ExtendStruct *ex ;

  if (isCheck && dd < 0) die ("dd %d < 0 in extendUpdate", dd) ;

  /* find highest place in stack less than or equal to dd */
  for (i = arrayMax(extendStack) ; i-- ; )
    { ex = arrp(extendStack, i, ExtendStruct) ;
      if (dd >= ex->d) break ;
      n0 += ex->n0 ; n1 += ex->n1 ;
    }
  arrayMax(extendStack) = i+1 ;
  if (dd > ex->d)	/* make a new exStack element */
    { ex = arrayp(extendStack, i+1, ExtendStruct) ;
      ex->d = dd ; ex->n0 = n0 ; ex->n1 = n1 ;
      ex->sum0 = ex[-1].sum0 + dd*n0 ; ex->sum1 = ex[-1].sum1 + dd*n1 ;
    }
  else		/* add these things from longer dd to current element */
    { ex->n0 += n0 ; ex->sum0 += dd*n0 ; ex->n1 += n1 ; ex->sum1 += dd*n1 ; }
  /* now add the current element on */
  if (x) { ++ex->n1 ; ex->sum1 += dd ; } else { ++ex->n0 ; ex->sum0 += dd ; }
  return ex ;
}

static void extendPrepare (PbwtCursor *u, int k)
{
  int j ;
  ExtendStruct *ex ;

  if (extendScoreSize != u->M)
    { if (extendScore0) { free (extendScore0) ; free (extendScore1) ; }
      extendScore0 = mycalloc (u->M+1, long) ; extendScore1 = mycalloc (u->M+1, long) ;
      extendScoreSize = u->M ;
    }

  /* first run through samples j == 1..M and add scores from top to j */
  extendStack = arrayReCreate (extendStack, 512, ExtendStruct) ;
  ex = arrayp (extendStack, 0, ExtendStruct) ;
  ex->d = 0 ; ex->n0 = 0 ; ex->n1 = 0 ; ex->sum0 = 1 ; ex->sum1 = 1 ;
  extendScore0[0] = 0 ; extendScore1[0] = 0 ;
  for (j = 0 ; j++ < u->M ; )
    { extendScore0[j] = ex->sum0 ; extendScore1[j] = ex->sum1 ;
      int bit = (j < u->M || u->d[j-1] > u->d[j]) ? u->d[j-1] : u->d[j] ;
      if (u->y[j-1]) extendScore1[j] += bit ; else extendScore0[j] += bit ;
      if (j < u->M) ex = extendUpdate (k - u->d[j], u->y[j-1]) ;
    }

  /* then reverse j == M-1..0 and add on scores from bottom to j */
  extendStack = arrayReCreate (extendStack, 512, ExtendStruct) ;
  ex = arrayp (extendStack, 0, ExtendStruct) ;
  ex->d = 0 ; ex->n0 = 0 ; ex->n1 = 0 ; ex->sum0 = 1 ; ex->sum1 = 1 ;
  for (j = u->M ; j-- > 0 ; )
    { extendScore0[j] += ex->sum0 ; extendScore1[j] += ex->sum1 ;
      int bit = (j || u->d[j+1] > u->d[j]) ? u->d[j+1] : u->d[j] ;
      if (u->y[j]) extendScore1[j] += bit ; else extendScore0[j] += bit ;
      if (j) ex = extendUpdate (k - u->d[j], u->y[j]) ;
    }
}

/* To be correct I should calculate the upper score up until but not including the 
   preceding j-1, and the lower score back until but not including the succeeding j, 
   because I don't know here the match lengths to the flanking j-1,j.  I can know them 
   during the dynamic programming, because I know the backtrace.  But this means 
   calculating and storing the dplus, dminus values induced by the backtrace, which adds 
   considerable complexity.  Above I assume that they are the greater of d[j] and d[j-1]
   (for j-1) and d[j] and d[j+1] (for j), which is a lower bound, taking care not to
   query out of bounds at the top and bottom.
*/

static inline double extendScore (int j, int x)
{
  double p = extendScore1[j] / (double) (extendScore0[j] + extendScore1[j]) ;
  return (x ? p : 1-p) ;
}

static int phaseExtend (int x0, int x1, PbwtCursor *uRef, int j0, 
			PhaseCell *pcOld, PhaseCell *pcNew, int k) 
/* returns the number of new back pointers created */
/* change at XXX for sampling or YYY for Viterbi */
/* this code uses extend structure above via extendScore, but did not perform well */
/* so replace with new phaseExtend() below, I hope */
{
  PhaseCell *old = &pcOld[j0] ;
  PhaseCell *new = &pcNew[pbwtCursorMap (uRef, x0, j0)] ;
  int j1New = pbwtCursorMap (uRef, x1, old->j1) ;
  double sBit = old->s * extendScore (j0, x0) * extendScore (old->j1, x1) ;
  if (isCheck && k == 51 && (j0 == 1389 || pbwtCursorMap(uRef,x0,j0) == 1370))
    printf ("at k %d extend j0 %d j1 %d s %.3g with x %d %d to j0New %d j1New %d sbit %.3g\n",
	    k, j0, old->j1, old->s, x0, x1, pbwtCursorMap (uRef, x0, j0), j1New, sBit) ;
  if (!new->s)		/* this is the first j to map here */
    { new->s = sBit ; new->sBest = sBit ;
      new->j1 = j1New ;
      new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
      return 1 ;
    }
  /* XXX for sampling calculate chance to pick this j rather than
     previous here */
  if (sBit > new->sBest) /* this is better than previous j that mapped here */
    { new->s = sBit ; 	 /* YYY originally +=, for Viterbi change to = */
      new->sBest = sBit ;
      new->j1 = j1New ;
      if (traceBackCell(new->back)->back == old->back) /* be careful - edge case! */
	{ traceBackCell(new->back)->value = x0 ; return 0 ; }
      else
	{ traceBackPrune (new->back, k, j0) ;
	  new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
	  return 1 ;
	}
    }
  /* new->s += sBit ;	   /* include if not Viterbi */
  return 0 ;
}
#endif // EXTEND0

#ifdef EXTEND1	/* alternate version minimises number of recombinations */
/* in this version pc->s is -1 minus the number of recombinations */
static char *extendMethodText = "EXTEND1" ;

static void extendPrepare (PbwtCursor *u, int k) {} 

static inline int phaseExtend (int x0, int x1, PbwtCursor *uRef, int j0,
			       PhaseCell *pcOld, PhaseCell *pcNew, int k) 
/* returns the number of new back pointers created */
{
  PhaseCell *old = &pcOld[j0] ;
  int j0New = pbwtCursorMap (uRef, x0, j0) ;
  PhaseCell *new = &pcNew[j0New] ;
  if (new->s && new->s > old->s) return 0 ; /* quick check whether this is possible */
  PhaseCell temp ;
  temp.s = old->s ;
  temp.dplus0 = pbwtCursorMapDplus (uRef, x0, j0, old->dplus0) ;
  temp.dminus0 = pbwtCursorMapDminus (uRef, x0, j0, old->dminus0) ;
  if (temp.dplus0 > old->dmax0 && temp.dminus0 > old->dmax0) /* need a recombination */
    { --temp.s ; temp.dmax0 = k ; } else temp.dmax0 = old->dmax0 ;
  temp.j1 = pbwtCursorMap (uRef, x1, old->j1) ;
  temp.dplus1 = pbwtCursorMapDplus (uRef, x1, old->j1, old->dplus1) ;
  temp.dminus1 = pbwtCursorMapDminus (uRef, x1, old->j1, old->dminus1) ;
  if (temp.dplus1 > old->dmax1 && temp.dminus1 > old->dmax1) /* need a recombination */
    { --temp.s ; temp.dmax1 = k ; } else temp.dmax1 = old->dmax1 ;

  if (!new->s)		/* this is the first j0 to map here */
    { *new = temp ;
      new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
      return 1 ;
    }
  else if (temp.s > new->s)	/* this is better than previous j that mapped here */
    { TraceBackIndex oldNewBack = new->back ;
      *new = temp ; 	/* originally +=, for Viterbi change to = */
      if (traceBackCell(oldNewBack)->back == old->back) /* be careful - edge case! */
	{ new->back = oldNewBack ; traceBackCell(new->back)->value = x0 ; return 0 ; }
      else
	{ traceBackPrune (oldNewBack, k, j0) ;
	  new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
	  return 1 ;
	}
    }
  else
    return 0 ;
}
#endif // EXTEND1

#ifdef EXTEND2	/* alternate version of EXTEND1 that minimises number of recombinations */
/* in this version we track the prefix interval (jmin,jmax) for matches to j0, j1 */
static char *extendMethodText = "EXTEND2" ;

static void extendPrepare (PbwtCursor *u, int k) {} 

static inline int phaseExtend (int x0, int x1, PbwtCursor *uRef, int j0,
			       PhaseCell *pcOld, PhaseCell *pcNew, int k) 
/* returns the number of new back pointers created */
{
  PhaseCell *old = &pcOld[j0] ;
  int j0New = pbwtCursorMap (uRef, x0, j0) ;
  PhaseCell *new = &pcNew[j0New] ;
  if (new->s && new->s > old->s) return 0 ; /* quick check whether this is possible */
  PhaseCell temp ;
  temp.s = old->s ;
  temp.j0min = pbwtCursorMap (uRef, x0, old->j0min) ;
  temp.j0max = pbwtCursorMap (uRef, x0, old->j0max) ;
  if (temp.j0min == temp.j0max) /* need a recombination */
    { --temp.s ; temp.j0min = x0 ? uRef->c-1 : 0 ; temp.j0max = x0 ? uRef->M : uRef->c ; }
  temp.j1 = pbwtCursorMap (uRef, x1, old->j1) ;
  temp.j1min = pbwtCursorMap (uRef, x1, old->j1min) ;
  temp.j1max = pbwtCursorMap (uRef, x1, old->j1max) ;
  if (temp.j1min == temp.j1max) /* need a recombination */
    { --temp.s ; temp.j1min = x1 ? uRef->c-1 : 0 ; temp.j1max = x1 ? uRef->M : uRef->c ; }

  if (!new->s)		/* this is the first j0 to map here */
    { *new = temp ;
      new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
      return 1 ;
    }
  else if (temp.s > new->s)	/* this is better than previous j that mapped here */
    { TraceBackIndex oldNewBack = new->back ;
      *new = temp ; 	/* originally +=, for Viterbi change to = */
      if (traceBackCell(oldNewBack)->back == old->back) /* be careful - edge case! */
	{ new->back = oldNewBack ; traceBackCell(new->back)->value = x0 ; return 0 ; }
      else
	{ traceBackPrune (oldNewBack, k, j0) ;
	  new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
	  return 1 ;
	}
    }
  else
    return 0 ;
}
#endif // EXTEND2

#ifdef EXTEND3	/* like EXTEND1,2 but maximises sum of lengths of maximal matches*/
/* in this version pc->s is 1 + sum of lengths of maximal matches */
static char *extendMethodText = "EXTEND3" ;

static void extendPrepare (PbwtCursor *u, int k) {} 

static inline int phaseExtend (int x0, int x1, PbwtCursor *uRef, int j0,
			       PhaseCell *pcOld, PhaseCell *pcNew, int k) 
/* returns the number of new back pointers created */
{
  PhaseCell *old = &pcOld[j0] ;
  int j0New = pbwtCursorMap (uRef, x0, j0) ;
  PhaseCell *new = &pcNew[j0New] ;
  PhaseCell temp ;	/* build the potential new in this  */
  temp.s = old->s ;
  temp.dplus0 = pbwtCursorMapDplus (uRef, x0, j0, old->dplus0) ;
  temp.dminus0 = pbwtCursorMapDminus (uRef, x0, j0, old->dminus0) ;
  int j ;
  if (old->dplus0 < old->dminus0 && old->dplus0 < temp.dplus0)
    for (j = j0 ; j < uRef->M ; )
      { temp.s += (k - old->dplus0) ;
	if (uRef->d[++j] > old->dplus0) break ;
      }
  else if (old->dminus0 < old->dplus0 && old->dminus0 < temp.dminus0)
    for (j = j0 ; j < uRef->M ; )
      { temp.s += (k - old->dminus0) ;
	if (uRef->d[--j] > old->dminus0) break ;
      }

  temp.j1 = pbwtCursorMap (uRef, x1, old->j1) ;
  temp.dplus1 = pbwtCursorMapDplus (uRef, x1, old->j1, old->dplus1) ;
  temp.dminus1 = pbwtCursorMapDminus (uRef, x1, old->j1, old->dminus1) ;
  if (old->dplus1 < old->dminus1 && old->dplus1 < temp.dplus1)
    for (j = j0 ; j < uRef->M ; )
      { temp.s += (k - old->dplus1) ;
	if (uRef->d[++j] > old->dplus1) break ;
      }
  else if (old->dminus1 < old->dplus1 && old->dminus1 < temp.dminus1)
    for (j = j0 ; j < uRef->M ; )
      { temp.s += (k - old->dminus1) ;
	if (uRef->d[--j] > old->dminus1) break ;
      }

  if (!new->s)		/* this is the first j0 to map here */
    { *new = temp ;
      new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
      return 1 ;
    }
  else if (temp.s > new->s)	/* this is better than previous j that mapped here */
    { TraceBackIndex oldNewBack = new->back ;
      *new = temp ; 	/* originally +=, for Viterbi change to = */
      if (traceBackCell(oldNewBack)->back == old->back) /* be careful - edge case! */
	{ new->back = oldNewBack ; traceBackCell(new->back)->value = x0 ; return 0 ; }
      else
	{ traceBackPrune (oldNewBack, k, j0) ;
	  new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
	  return 1 ;
	}
    }
  else
    return 0 ;
}
#endif // EXTEND3

#ifdef EXTEND4	/* based on diff from neighbours based on -log generative model */
/* in this version pc->s is -1 - sum of lengths of mismatched neighbours */
static char *extendMethodText = "EXTEND4" ;

static void extendPrepare (PbwtCursor *u, int k) {} 

static inline int phaseExtend (int x0, int x1, PbwtCursor *uRef, int j0,
			       PhaseCell *pcOld, PhaseCell *pcNew, int k) 
/* returns the number of new back pointers created */
{
  PhaseCell *old = &pcOld[j0] ;
  int j0New = pbwtCursorMap (uRef, x0, j0) ;
  PhaseCell *new = &pcNew[j0New] ;
  PhaseCell temp ;	/* build the potential new in this  */
  temp.s = old->s ;
  temp.dplus0 = pbwtCursorMapDplus (uRef, x0, j0, old->dplus0) ;
  temp.dminus0 = pbwtCursorMapDminus (uRef, x0, j0, old->dminus0) ;
  float dS = 0 ;
  if (j0) 
    { if (x0 == uRef->y[j0-1]) dS += (k-old->dminus0) ; else dS -= (k-old->dminus0) ; }
  if (j0 < uRef->M) 
    { if (x0 == uRef->y[j0]) dS += (k-old->dplus0) ; else dS -= (k-old->dplus0) ; }
  if (dS < 0) temp.s += dS ;

  temp.j1 = pbwtCursorMap (uRef, x1, old->j1) ;
  temp.dplus1 = pbwtCursorMapDplus (uRef, x1, old->j1, old->dplus1) ;
  temp.dminus1 = pbwtCursorMapDminus (uRef, x1, old->j1, old->dminus1) ;
  dS = 0 ;
  if (old->j1) 
    { if (x1 == uRef->y[old->j1-1]) dS += (k-old->dminus1) ; else dS -= (k-old->dminus1) ; }
  if (old->j1 < uRef->M) 
    { if (x1 == uRef->y[old->j1]) dS += (k-old->dplus1) ; else dS -= (k-old->dplus1) ; }
  if (dS < 0) temp.s += dS ;

  if (!new->s)		/* this is the first j0 to map here */
    { *new = temp ;
      new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
      return 1 ;
    }
  else if (temp.s > new->s)	/* this is better than previous j that mapped here */
    { TraceBackIndex oldNewBack = new->back ;
      *new = temp ; 	/* originally +=, for Viterbi change to = */
      if (traceBackCell(oldNewBack)->back == old->back) /* be careful - edge case! */
	{ new->back = oldNewBack ; traceBackCell(new->back)->value = x0 ; return 0 ; }
      else
	{ traceBackPrune (oldNewBack, k, j0) ;
	  new->back = (x0 == x1) ? old->back : traceBackCreate (old->back, x0, j0, k) ;
	  return 1 ;
	}
    }
  else
    return 0 ;
}
#endif // EXTEND4

/********** finally the main function ********/

static PBWT *referencePhase4 (PBWT *pOld, PBWT *pRef)
{
  fprintf (logFile, "Reference phase with extension method %s\n", extendMethodText) ;
  int i, j, jq, k ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  uchar *xOld = myalloc (pOld->M, uchar) ;
  PbwtCursor *uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  traceBackInitialise() ;
  PhaseCell **pcOld = myalloc (pOld->M, PhaseCell*) ;
  PhaseCell **pcNew = myalloc (pOld->M, PhaseCell*) ;
  for (jq = 0 ; jq < pOld->M ; jq += 2)
    { pcOld[jq] = mycalloc (pRef->M+1, PhaseCell) ;
      pcOld[jq][0].back = 0 ;
#if defined(EXTEND0) || defined(EXTEND3)
      pcOld[jq][0].s = 1.0 ;
#endif
#if defined(EXTEND1) || defined(EXTEND2) || defined(EXTEND4)
      pcOld[jq][0].s = -1.0 ;
#endif
      pcNew[jq] = myalloc (pRef->M+1, PhaseCell) ;
    }
  /* some declarations for checks */
  int *checkLiveSum = mycalloc (pOld->M, int) ; 
  double *checkLikeSum = mycalloc (pOld->M, double) ; 
  int *checkHets = mycalloc (pOld->M, int) ;

  /* here is the main loop through the sites */
  for (k = 0 ; k < pOld->N ; ++k)
    { for (j = 0 ; j < pOld->M ; ++j) xOld[uOld->a[j]] = uOld->y[j] ;
#ifdef QUERY0_IS_REF0
      if (isCheck)
	{ PhaseCell *pc ; TraceBackCell *tc ;
	  int a0, a1 ;
	  for (j = 0 ; j < pRef->M ; ++j) 
	    { if (uRef->a[j] == 0) a0 = j ;
	      if (uRef->a[j] == 1) a1 = j ;
	    }
	  PhaseCell *pc0 = &pcOld[0][a0], *pc1 = &pcOld[0][a1] ;
	  printf ("at k %4d x0 %d x1 %d", k, xOld[0], xOld[1]) ;
	  printf ("\n  0 j %4d s %10.3g", a0, pc0->s) ; traceBackReport (pc0->back) ;
	  printf ("\n  1 j %4d s %10.3g", a1, pc1->s) ; traceBackReport (pc1->back) ;
	  for (j = 0, pc = pcOld[0] ; j < pRef->M ; ++j, ++pc) 
	    if (pc->s > pc0->s && pc->s > pc1->s) 
	      { printf ("\n    j %4d s %10.3g", j, pc->s) ; traceBackReport (pc->back) ; 
		die ("done") ;
	      }
	  putchar ('\n') ; 
	}
#endif
      if (isCheck && !(k % 100)) 
	{ printf ("check k %d", k) ; int nS = 0 ; PhaseCell *pc ;
	  for (j = 0, pc = pcOld[0] ; j < pRef->M ; ++j, ++pc) if (pc->s) ++nS ;
	  printf (" nS %d\n", nS) ;
	}
      extendPrepare (uRef, k) ;
      pbwtCursorCalculateU (uRef) ; /* before pbwtCursorMap()s in phaseExtend() */
      /* update each query sample in turn */
      for (jq = 0 ; jq < pOld->M ; jq +=2)
	{ bzero (pcNew[jq], (pRef->M+1)*sizeof(PhaseCell)) ; /* clear pcNew */
	  int x0 = xOld[jq], x1 = xOld[jq+1] ;
	  if (x0 != x1) checkHets[jq]++ ;
	  for (j = 0 ; j <= pRef->M ; ++j)
	    if (pcOld[jq][j].s)
	      { int countBack ;
		if (x0 == x1)
		  countBack = phaseExtend (x0, x1, uRef, j, pcOld[jq], pcNew[jq], k) ;
		else 
		  countBack = phaseExtend (x0, x1, uRef, j, pcOld[jq], pcNew[jq], k) + 
		              phaseExtend (x1, x0, uRef, j, pcOld[jq], pcNew[jq], k) ;
		/* now resolve traceback counts on pcOld[jq][j].back */
		if (countBack == 0) 
		  traceBackPrune (pcOld[jq][j].back, k, j) ;
		else if (countBack == 2 && pcOld[jq][j].back) /* don't increment tb cell 0 */
		  traceBackCell(pcOld[jq][j].back)->count++ ;
		/* else leave count at default of 1 */
		++checkLiveSum[jq] ;
	      }
	  /* now update pcOld from pcNew to be ready for next column */
 	  double sum = 0 ; for (j = 0 ; j <= pRef->M ; ++j) sum += pcNew[jq][j].s ;
	  if (!sum) die ("sum is 0 at k %d jq %d", k, jq) ;
	  checkLikeSum[jq] += log (sum) ;
	  for (j = 0 ; j <= pRef->M ; ++j)
	    { pcOld[jq][j] = pcNew[jq][j] ;
#ifdef EXTEND0
	      if (!isCheck) pcOld[jq][j].s /= sum ; /* normalise */
#endif
	    }
	}
#ifdef TRACEBACK_DEBUG
      if (isCheck) traceBackCheck (jq, k, pcOld[0], pRef->M) ; /* check traceback counters */
      printf ("%d%d    %d", xOld[0], xOld[1], uRef->y[0]) ;
      for (j = 1 ; j <= pRef->M ; ++j) printf ("      %d", uRef->y[j]) ; putchar ('\n') ;
      for (j = 0 ; j <= pRef->M ; ++j) printf (" %6.1f", pcOld[0][j].s) ; putchar ('\n') ;
      for (j = 0 ; j <= pRef->M ; ++j) printf (" %6d", pcOld[0][j].j1) ; putchar ('\n') ;
#endif

      pbwtCursorForwardsRead (uOld) ;
      pbwtCursorForwardsReadAD (uRef, k) ;
    }

  fprintf (logFile, "traceBackHeap final %ld, max %ld\n", 
	   arrayMax(traceBackHeap)-arrayMax(traceBackFreeStack), arrayMax(traceBackHeap)) ;

  /* now do the traceback - we write first into the reverse pbwt for pNew */
  /* first find highest score in last column, to start at */
  TraceBackIndex *tb = myalloc (pOld->M, TraceBackIndex) ;
  for (jq = 0 ; jq < pOld->M ; jq += 2)
    { double sMax = -1e20 ; int jMax ; 
      for (j = 0 ; j <= pRef->M ; ++j) 
	if (pcOld[jq][j].s && pcOld[jq][j].s > sMax) { sMax = pcOld[jq][j].s ; jMax = j ; }
      tb[jq] = pcOld[jq][jMax].back ;
    }
  /* now trace back from there, building pNew backwards */
  PBWT *pNew = pbwtCreate (pOld->M, pOld->N) ; /* need to initialist aRstart */
  pNew->aRstart = myalloc (pNew->M, int) ; for (j=0 ; j < pNew->M ; ++j) pNew->aRstart[j] = j ;
  PbwtCursor *uNewR = pbwtCursorCreate (pNew, FALSE, TRUE) ;
  uchar *xNew = myalloc (pNew->M, uchar) ;
  for (k = pOld->N ; k-- ; )
    { pbwtCursorReadBackwards (uOld) ;
      for (j = 0 ; j < pOld->M ; ++j) xOld[uOld->a[j]] = uOld->y[j] ;
      for (jq = 0 ; jq < pOld->M ; jq += 2)
	if (xOld[jq] == xOld[jq+1])
	  { xNew[jq] = xOld[jq] ; xNew[jq+1] = xOld[jq+1] ; }
	else
	  { TraceBackCell *tc = traceBackCell(tb[jq]) ;
	    xNew[jq] = tc->value ; xNew[jq+1] = 1 - xNew[jq] ;
#ifdef TRACEBACK_DEBUG
	    if (k != tc->k) die ("k mismatch during traceback k %d tc->k %d jq %d") ;
	    if (isCheck) 
	      printf ("traceBack k %d jq %d tb %d tc->value %d ->k %d ->j %d ->back %d\n",
		      k, jq, tb[jq], tc->value, tc->k, tc->j, tc->back) ;
#endif
	    if (!tb[jq]) die ("premature end of trace back at k %d, jq %d", k, jq) ;
	    tb[jq] = tc->back ;
	  }
      for (j = 0 ; j < pNew->M ; ++j) uNewR->y[j] = xNew[uNewR->a[j]] ;
      pbwtCursorWriteForwards (uNewR) ;
    }
  for (jq = 0 ; jq < pOld->M ; jq += 2) if (tb[jq]) die ("trace back incomplete jq %d", jq) ;
  /* copy final sort order into pNew->aRend and ->aFstart */
  pNew->aRend = myalloc (pNew->M, int) ; memcpy (pNew->aRend, uNewR->a, pNew->M*sizeof(int)) ;
  memcpy (pNew->aFstart, uNewR->a, pNew->M * sizeof(int)) ;
  /* finally make the forwards pbwt for pNew */
  PbwtCursor *uNewF = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  for (k = 0 ; k < pOld->N ; k++)
    { pbwtCursorReadBackwards (uNewR) ;
      for (j = 0 ; j < pNew->M ; ++j) xOld[uNewR->a[j]] = uNewR->y[j] ;
      for (j = 0 ; j < pNew->M ; ++j) uNewF->y[j] = xOld[uNewF->a[j]] ;
      pbwtCursorWriteForwards (uNewF) ;
    }
  pNew->aFend = myalloc (pNew->M, int) ; memcpy (pNew->aFend, uNewF->a, pNew->M*sizeof(int)) ;

  /* reporting */
  if (isCheck)
    for (jq = 0 ; jq < pOld->M ; jq +=2 )
      fprintf (logFile, "jq %d, nHets %d, liveAv %.2f, likeAv %.2f\n", jq, 
	       checkHets[jq], checkLiveSum[jq]/(double)pRef->N, exp(checkLikeSum[jq]/pRef->N)) ;

  /* cleanup */
  pbwtCursorDestroy (uRef) ; pbwtCursorDestroy (uOld) ;
  pbwtCursorDestroy (uNewF) ; pbwtCursorDestroy (uNewR) ;
  for (jq = 0 ; jq < pOld->M ; jq +=2)
    { free (pcOld[jq]) ; free (pcNew[jq]) ; }
  free (pcOld) ; free (pcNew) ;
  free (xOld) ; free (xNew) ;
  free (checkLiveSum) ; free (checkLikeSum) ; free (checkHets) ;

  return pNew ;
}

/************* top level selector for referencePhase *****************/

PBWT *referencePhase (PBWT *pOld, char *fileNameRoot)
{
  fprintf (logFile, "phase against reference %s\n", fileNameRoot) ;
  if (pOld->M % 2) die ("phase requires that M = %d is even", pOld->M) ;
  if (!pOld || !pOld->yz || !pOld->sites) 
    die ("referencePhase called without existing pbwt with sites") ;

  PBWT *pRef = pbwtReadAll (fileNameRoot) ;
  if (!pRef->sites) die ("new pbwt %s in referencePhase has no sites", fileNameRoot) ;
  if (strcmp(pOld->chrom,pRef->chrom))
    die ("mismatching chrom in referencePhase: old %s, ref %s", pOld->chrom, pRef->chrom) ;

  /* reduce both down to the intersecting sites */
  pOld = pbwtSelectSites (pOld, pRef->sites, FALSE) ;
  pRef = pbwtSelectSites (pRef, pOld->sites, FALSE) ;
  if (!pOld->N) die ("no overlapping sites in referencePhase") ;

  fprintf (logFile, "Phase preliminaries: ") ; timeUpdate(logFile) ;
  PBWT *pNew = referencePhase4 (pOld, pRef) ;
  fprintf (logFile, "Phasing complete: ") ; timeUpdate(logFile) ;
  fprintf (logFile, "After phasing: ") ; phaseCompare (pNew, pOld) ;

  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->sites = pOld->sites ; pOld->sites = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtDestroy (pRef) ;
  return pNew ;
}

/********** impute using maximal matches - use matchSequencesSweep() to find them *********/

static double **pImp = 0 ;	/* use with stats to store p values */

int matchSegmentCompare (const void *a, const void *b) /* for sorting */
{
  return ((MatchSegment*)a)->start - ((MatchSegment*)b)->start ;
}

/* overload a high-end bit of the ->end field of MatchSegment to record if match is sparse */
static int SPARSE_BIT = 1 << 30 ;
static int SPARSE_MASK = (1 << 30) - 1 ;

static void reportMatchSparse (int iq, int jRef, int start, int end, BOOL isSparse)
{
  MatchSegment *ms = arrayp (maxMatch[iq], arrayMax(maxMatch[iq]), MatchSegment) ;
  ms->jRef = jRef ; ms->start = start ; ms->end = end ;
  if (isSparse) ms->end |= SPARSE_BIT ;
}

static PBWT *referenceImpute3 (PBWT *pOld, PBWT *pRef, PBWT *pFrame, 
			       int nSparse, double fSparse)
/* Require pOld and pFrame to have the same sites, a subset of sites of pRef, */
/* and pRef and pFrame to have the same samples. */
/* Impute only missing sites in pRef if pOld == pFrame. */
/* Added nSparse to allow also matching at sparse positions */
{
  int i, j, k ;

  fprintf (logFile, "Reference impute using maximal matches: ") ;
  if (nSparse > 1) fprintf (logFile, "(nSparse = %d, fSparse = %.2f) ", nSparse, fSparse) ;

  /* build the array of maximal matches into pFrame for each sequence in pOld */
  maxMatch = myalloc (pOld->M, Array) ;
  for (j = 0 ; j < pOld->M ; ++j) maxMatch[j] = arrayCreate (1024, MatchSegment) ;
  if (pOld == pFrame)		/* self-imputing - no sparse option yet */
    matchMaximalWithin (pFrame, reportMatch) ;
  else
    matchSequencesSweep (pFrame, pOld, reportMatch) ;
    //matchSequencesSweepSparse (pFrame, pOld, nSparse, reportMatchSparse) ;

  for (j = 0 ; j < pOld->M ; ++j)	/* add terminating element to arrays */
    { if (nSparse > 1) /* can't guarantee order of sparse segments */
#ifdef MERGESORT	
	if (mergesort (arrp(maxMatch[j], 0, MatchSegment), arrayMax(maxMatch[j]), 
		       sizeof(MatchSegment), matchSegmentCompare))
	  die ("error %d in mergesort", errno) ;
      /* mergesort() because they are close to being already sorted */
#else
        qsort (arrp(maxMatch[j], 0, MatchSegment), arrayMax(maxMatch[j]), 
	       sizeof(MatchSegment), matchSegmentCompare) ;
#endif
      if (isCheck) fprintf (logFile, "%ld matches found to query %d\n", 
			    arrayMax(maxMatch[j]), j) ;
      /* add an end marker */
      MatchSegment *ms = arrayp(maxMatch[j],arrayMax(maxMatch[j]),MatchSegment) ;
      ms->jRef = ms[-1].jRef ; ms->end = pOld->N+1 ; ms->start = pOld->N ;
    }

  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uRef = pbwtCursorCreate (pRef, TRUE, TRUE) ;
  PBWT *pNew = pbwtCreate (pOld->M, pRef->N) ;	/* this will hold the imputed sequence */
  pNew->isRefFreq = TRUE ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  uchar *x = myalloc (pOld->M, uchar) ;     /* uNew->y values in original sort order */
  // unused: double *p = myalloc (pOld->M, double) ;   /* estimated prob of uNew->y in original sort order */
  int *aRefInv = myalloc (pRef->M, int) ;   /* holds the inverse mapping from uRef->a[i] -> i */
  int *firstSeg = mycalloc (pOld->M, int) ; /* position in maxMatch to start looking at */
  int nConflicts = 0 ;
  uchar *missing = mycalloc (pOld->M, uchar) ;
  double *xDosage = myalloc (pOld->M, double), *yDosage = myalloc (pOld->M, double) ;

  pNew->dosageOffset = arrayReCreate (pNew->dosageOffset, pRef->N, long) ;
  pNew->zDosage = arrayReCreate (pNew->zDosage, pRef->N*16, uchar) ;

  int kOld = 0, kRef = 0 ;
  while (kRef < pRef->N)
    { if (arrp(pRef->sites,kRef,Site)->x == arrp(pFrame->sites,kOld,Site)->x
	  && arrp(pRef->sites,kRef,Site)->varD == arrp(pFrame->sites,kOld,Site)->varD)
	{ 
      if ( pOld != pFrame )
      {
          if (!pOld->missingOffset || !arr(pOld->missingOffset, kOld, long)) bzero (missing, pOld->M) ;
          else unpack3 (arrp(pOld->zMissing,arr(pOld->missingOffset,kOld,long), uchar), pOld->M, missing, 0) ;
      }
      pbwtCursorForwardsRead (uOld) ; ++kOld ;
	  for (j = 0 ; j < pOld->M ; ++j)
	    while (kOld > (arrp(maxMatch[j],firstSeg[j],MatchSegment)->end & SPARSE_MASK)) ++firstSeg[j] ;
	}
      else
        arrp(pRef->sites,kRef,Site)->isImputed = TRUE ;
      for (i = 0 ; i < pRef->M ; ++i) aRefInv[uRef->a[i]] = i ;
      double psum = 0, xsum = 0, pxsum = 0 ; int n = 0 ;
      arrp(pRef->sites,kRef,Site)->refFreq = (uRef->M - uRef->c) / (double) pRef->M ;
      if (pOld == pFrame)	/* find which samples are missing at this site */
	{ if (!arr(pRef->missingOffset, kRef, long)) bzero (missing, pRef->M) ;
	  else unpack3 (arrp(pRef->zMissing,arr(pRef->missingOffset,kRef,long), uchar), 
			pRef->M, missing, 0) ;
	}
      for (j = 0 ; j < pOld->M ; ++j)
	{ if (pOld == pFrame && !missing[j]) /* don't impute - copy from ref */
	    { x[j] = uRef->y[aRefInv[j]] ;
	      continue ;
	    }
	  /* impute from overlapping matches */
	  double bit = 0, sum = 0, score = 0 ;
	  MatchSegment *m = arrp(maxMatch[j],firstSeg[j],MatchSegment) ;
	  MatchSegment *mStop = arrp(maxMatch[j],arrayMax(maxMatch[j]),MatchSegment) ;
	  while (m->start <= kOld && m->end >= kOld && m < mStop)
	    { bit = (kOld - m->start) * ((m->end & SPARSE_MASK) - kOld + 1) ;
	      if (m->end & SPARSE_BIT) bit *= fSparse ;
	      if (bit > 0)
		{ sum += bit ;
		  if (uRef->y[aRefInv[m->jRef]]) score += bit ;
		}
	      ++m ;
	    }
	  if (sum == 0) 
	    { x[j] = arrp(pRef->sites,kRef,Site)->refFreq > 0.5 ? 1 : 0 ;
	      xDosage[j] = arrp(pRef->sites,kRef,Site)->refFreq ;
	      if (isStats) pImp[kRef][j] = xDosage[j] ;
	      ++nConflicts ;
	    }
	  else 
	    { double p = score/sum ;
	      x[j] = (p > 0.5) ? 1 : 0 ;
	      xDosage[j] = p ;
	      psum += p ;
	      xsum += x[j] ;
	      pxsum += p*x[j] ;
	      if (isStats) pImp[kRef][j] = p ;
	      ++n ;
	    }
	}
	  
      for (j = 0 ; j < pOld->M ; ++j) uNew->y[j] = x[uNew->a[j]] ; /* transfer to uNew */
      /* need to sort the dosages into uNew cursor order as well */
      for (j = 0 ; j < pOld->M ; ++j) 
      {
        /* make sure the dosages of the typed markers are {0,1} */
        if ( arrp(pRef->sites,kRef,Site)->isImputed!=TRUE && !missing[uNew->a[j]] ) yDosage[j] = x[uNew->a[j]] ? 1 : 0;
        else yDosage[j] = xDosage[uNew->a[j]] ;
      }
      pbwtDosageStore (pNew, yDosage, kRef) ;
      pbwtCursorWriteForwards (uNew) ;
      
      if (n) 
	{ psum /= n ; xsum /= n ; pxsum /= n ;
	  double varianceProduct = psum*(1.0-psum)*xsum*(1.0-xsum) ;
	  if (varianceProduct)
	    arrp(pRef->sites,kRef,Site)->imputeInfo = 
	      (pxsum - psum*psum) / sqrt (varianceProduct) ;
	  else
	    arrp(pRef->sites,kRef,Site)->imputeInfo = 1.0 ;
	}
      pbwtCursorForwardsRead (uRef) ; ++kRef ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  if (nConflicts) fprintf (logFile, "%d times where no overlapping matches because query does not match any reference - set imputed value to 0\n", nConflicts) ;

  pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uRef) ; pbwtCursorDestroy (uNew) ;
  free (aRefInv) ; free (firstSeg) ;
  for (j = 0 ; j < pOld->M ; ++j) arrayDestroy (maxMatch[j]) ; free (maxMatch) ;
  free (xDosage) ; free (yDosage) ; if (missing) free (missing) ; free (x) ;
  return pNew ;
}

/*********************************************************************/

  PBWT *referenceImpute (PBWT *pOld, char *fileNameRoot, int nSparse, double fSparse)
{
  if ( !pOld->N ) return pOld;  /* empty file, don't segfault */

  /* Preliminaries */
  fprintf (logFile, "impute against reference %s\n", fileNameRoot) ;
  if (!pOld || !pOld->yz || !pOld->sites) 
    die ("referenceImpute called without existing pbwt with sites") ;
  PBWT *pRef = pbwtReadAll (fileNameRoot) ;
  if (!pRef->sites) die ("new pbwt %s in referencePhase has no sites", fileNameRoot) ;
  if (strcmp(pOld->chrom,pRef->chrom))
    die ("mismatching chrom in referenceImpute: old %s, new %s", pRef->chrom, pOld->chrom) ;

  /* identify the intersecting sites */
  PBWT *pFrame = pbwtSelectSites (pRef, pOld->sites, TRUE) ; /* keep the full ref to impute to */
  pFrame->isX = pRef->isX ; pFrame->isY = pRef->isY ;
  if (pFrame->N == pRef->N)
    { fprintf (logFile, "No additional sites to impute in referenceImpute\n") ;
      pbwtDestroy (pFrame) ; pbwtDestroy (pRef) ;
      return pOld ;
    }
  pbwtBuildReverse (pFrame) ;	/* we need the reverse reference pbwt below */
  pOld = pbwtSelectSites (pOld, pRef->sites, FALSE) ;
  if (!pOld->N) die ("no overlapping sites in referenceImpute") ;
  if (!pOld->aFend) die ("pOld has no aFend in referenceImpute - your pbwt was made by a previous version of the code; buildReverse and resave the forwards pbwt") ;

  fprintf (logFile, "Imputation preliminaries: ") ; timeUpdate(logFile) ;

  if (isStats)
    { pImp = myalloc (pRef->N, double*) ;
      int k ; for (k = 0 ; k < pRef->N ; ++k) pImp[k] = myalloc (pOld->M, double) ;
    }

  PBWT *pNew = referenceImpute3 (pOld, pRef, pFrame, nSparse, fSparse) ;
  pNew->sites = pRef->sites ; pRef->sites = 0 ; 
  pNew->chrom = pRef->chrom ; pRef->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;

  if (isStats)
    { int k, j, ff ; long his[20][10] ;
      bzero (his, 20*10*sizeof(long)) ;
      for (k = 0 ; k < pRef->N ; ++k)
	{ double f = arrp(pNew->sites,k,Site)->refFreq ;
	  for (ff = 0 ; f*100 > fBound[ff] ; ++ff) ;
	  for (j = 0 ; j < pNew->M ; ++j) ++his[ff][(int)(pImp[k][j]*10)] ;
	}
      for (ff = 0 ; ff < 17 ; ++ff)
	{ fprintf (logFile, "%5.1f", fBound[ff]) ;
	  double tot = 0.0 ; for (j = 10 ; j-- ;)  tot += his[ff][j] ;
	  for (j = 0 ; j < 10 ; ++j) 
	    fprintf (logFile, " %8.5f", his[ff][j]/tot) ;
	  fprintf (logFile, "\n") ;
	}
    }

  pbwtDestroy (pOld) ; pbwtDestroy (pFrame) ; pbwtDestroy (pRef) ;
  return pNew ;
}

/*********************************************************************/

PBWT *imputeMissing (PBWT *pOld)
/* current strategy is for HRC: use framework of sites for which we have complete data */
{
  int k ;

  if (!pOld->missingOffset)
    { warn ("imputeMissing called but can't find missing data\n") ;
      return pOld ;
    }

  if (isStats)
    { int *nMiss = mycalloc (10,int), n0, i ;
      uchar *miss = myalloc (pOld->M, uchar) ;
      for (k = 0 ; k < pOld->N ; ++k)
	if (arr(pOld->missingOffset,k,long)) 
	  { unpack3 (arrp(pOld->zMissing,arr(pOld->missingOffset,k,int), uchar), 
		     pOld->M, miss, &n0) ;
	    n0 = pOld->M - n0 ;
	    if (isCheck)
	      { printf ("missing at %d: offset %ld, n0 %d", 
			k, arr(pOld->missingOffset,k,long), n0) ;
		putchar ('\n') ;
	      }
	    for (i = 0 ; n0 > 0 ; ++i) { ++nMiss[i] ; n0 /= 10 ; }
	  }
      n0 = 1 ; 
      for (i = 0 ; nMiss[i] ; ++i)
	{ printf ("sites with missing >= %d: %d\n", n0, nMiss[i]) ;
	  n0 *= 10 ;
	}
      free (miss) ; free (nMiss) ;
    }

  /* first build frame */
  Array completeSites = arrayCreate (pOld->M, Site) ;
  for (k = 0 ; k < pOld->N ; ++k) 
    if (!arr(pOld->missingOffset,k,long)) 
      array(completeSites,arrayMax(completeSites),Site) = arr(pOld->sites,k,Site) ;
  PBWT *pFrame = pbwtSelectSites (pOld, completeSites, TRUE) ;
  arrayDestroy (completeSites) ;

  /* then impute, using a special mode of impute2() */
  PBWT *pNew = referenceImpute3 (pFrame, pOld, pFrame, 1, 0) ;
  pNew->sites = pOld->sites ; pOld->sites = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;
  pbwtDestroy (pOld) ; pbwtDestroy (pFrame) ;
  return pNew ;
}

/*********************************************************************/

void genotypeCompare (PBWT *p, char *fileNameRoot)
{
  fprintf (logFile, "compare genotypes to reference %s\n", fileNameRoot) ;
  if (!p || !p->yz || !p->sites) 
    die ("genotypeCompare called without existing pbwt with sites") ;
  PBWT *pRef = pbwtReadAll (fileNameRoot) ;
  if (strcmp(p->chrom,pRef->chrom)) die ("mismatch chrom %s to ref %f", p->chrom, pRef->chrom) ;
  if (!pRef->sites) die ("new pbwt %s in genotypeCompare has no sites", fileNameRoot) ;
  if (p->M != pRef->M) die ("mismatch of old M %d to ref M %d", p->M, pRef->M) ;
  if (p->N == pRef->N)
      genotypeComparePbwt (p, pRef) ;
  else		      /* reduce both down to the intersecting sites */
    { warn ("mismatch of old N %d to ref N %d", p->N, pRef->N) ;
      PBWT *pFrame = pbwtSelectSites (p, pRef->sites, TRUE) ;
      pRef = pbwtSelectSites (pRef, p->sites, FALSE) ;
      if (!pFrame->N) die ("no overlapping sites in genotypeCompare") ;
      genotypeComparePbwt (pFrame, pRef) ;
      pbwtDestroy(pFrame) ;
    }
  
  pbwtDestroy (pRef) ;
}

static void genotypeComparePbwt (PBWT *p, PBWT *q)
{
  int i, j, k, ff ;
  long n[17][9] ; for (i = 17 ; i-- ;) for (j = 9 ; j-- ;) n[i][j] = 0 ;
  double fsum[17] ; for (i = 17 ; i-- ;) fsum[i] = 0.0 ;
  int nsum[17] ; for (i = 17 ; i-- ;) nsum[i] = 0 ;
  double isum[17] ; for (i = 17 ; i-- ;) isum[i] = 0.0 ;
  int ni[17] ; for (i = 17 ; i-- ;) ni[i] = 0 ;
  long *ns = mycalloc (9*p->M, long) ;
  BOOL isRefFreq = FALSE ;
  BOOL isDosage = p->dosageOffset ? TRUE : FALSE ;
  double *dos = 0 ;
  long nd[12] ; for (i = 12 ; i-- ;) nd[i] = 0 ;
  long nd1[12] ; for (i = 12 ; i-- ;) nd1[i] = 0 ;

  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  PbwtCursor *uq = pbwtCursorCreate (q, TRUE, TRUE) ;
  uchar *xp = myalloc (p->M, uchar), *xq = myalloc (q->M, uchar) ;
  for (k = 0 ; k < p->N ; ++k)
    { double f = (p->M - up->c) / (double)p->M ; 
      if (arrp (p->sites,k,Site)->refFreq) { f = arrp (p->sites,k,Site)->refFreq ; isRefFreq = TRUE ; }
      for (ff = 0 ; f*100 > fBound[ff] ; ) ++ff ; fsum[ff] += f*100 ; ++nsum[ff] ;
      if (arrp(p->sites,k,Site)->imputeInfo < 1.0) { isum[ff] += arrp(p->sites,k,Site)->imputeInfo ; ++ni[ff] ; }
      for (j = 0 ; j < p->M ; ++j) { xp[up->a[j]] = up->y[j] ; xq[uq->a[j]] = uq->y[j] ; }
      if (isDosage) dos = pbwtDosageRetrieve (p, up, dos, k) ;
      for (j = 0 ; j < p->M ; j += 2) 
	{ i = 3*(xp[j]+xp[j+1]) + xq[j]+xq[j+1] ;
	  ++n[ff][i] ;
	  ++ns[9*j + i] ;
	  if (isDosage) 
	    { int id ; 
	      if (dos[j] == 0.0) id = 0 ; else if (dos[j] == 1.0) id = 11 ; else id = 1 + (int)(dos[j]*10.0) ;
	      ++nd[id] ; if (xp[j]) ++nd1[id] ;
	      if (dos[j+1] == 0.0) id = 0 ; else if (dos[j+1] == 1.0) id = 11 ; else id = 1 + (int)(dos[j+1]*10.0) ;
	      ++nd[id] ; if (xp[j+1]) ++nd1[id] ;
	    }
	}
      pbwtCursorForwardsRead (up) ; pbwtCursorForwardsRead (uq) ;
    }
  pbwtCursorDestroy (up) ; pbwtCursorDestroy (uq) ;

  /* report */

  if (isRefFreq) printf ("Genotype comparison results split on reference frequencies\n") ;
  else printf ("Genotype comparison results split on sample frequencies\n") ;
  for (ff = 0 ; ff < 17 ; ++ff)
    { double tot = 0, xbar, ybar, r2, x2, y2 ;
      printf ("%-5.1f\t%-7.3f", fBound[ff], nsum[ff] ? fsum[ff]/nsum[ff] : 0.0) ;
      for (i = 0 ; i < 9 ; ++i) { printf ("\t%ld ", n[ff][i]) ; tot += n[ff][i] ; }
      if (tot) 
	{ xbar = (n[ff][3] + n[ff][4] + n[ff][5] + 2*(n[ff][6] + n[ff][7] + n[ff][8])) / tot ;
	  x2 = (n[ff][3] + n[ff][4] + n[ff][5] + 4*(n[ff][6] + n[ff][7] + n[ff][8])) / tot ;
	  ybar = (n[ff][1] + n[ff][4] + n[ff][7] + 2*(n[ff][2] + n[ff][5] + n[ff][8])) / tot ;
	  y2 = (n[ff][1] + n[ff][4] + n[ff][7] + 4*(n[ff][2] + n[ff][5] + n[ff][8])) / tot ;
	  r2 = (n[ff][4] + 2*(n[ff][5] + n[ff][7]) + 4*n[ff][8]) / tot ;
	  r2 = (r2 - xbar*ybar)/sqrt((x2 - xbar*xbar)*(y2 - ybar*ybar)) ;
	  printf ("\tx,y,r2\t%.4f\t%.4f\t%.4f", xbar, ybar, r2) ;
	  if (ni[ff]) printf ("\t info %.4f", isum[ff]/ni[ff]) ;
	}
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

  printf ("Genotype accuracy distribution across samples\n") ;
  if (hist[100]) printf ("%d samples with r2 == 1.0\n", hist[100]) ;
  for (i = 100 ; i-- ; ) 
    if (hist[i])
      printf ("%d samples with %.2f <= r2 < %.2f\n", hist[i], (i-1)*0.01, i*0.01) ;

  if (isDosage)
    { printf ("Dosage accuracy (currently at haplotype level)\n") ;
      printf ("0.00  %.3f  %ld\n", nd[0] ? nd1[0]/(double)nd[0] : 0.0, nd[0]) ;
      for (i = 1 ; i < 11 ; ++i)
	printf ("%.2f  %.3f  %ld\n", 0.1*(i-0.5), nd[i] ? nd1[i]/(double)nd[i] : 0.0, nd[i]) ;
      printf ("1.00  %.3f  %ld\n", nd[11] ? nd1[11]/(double)nd[11] : 0.0, nd[11]) ;
    }
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

  fprintf (logFile, "corruptSites with pSite %f, pChange %f changes %.4f of values\n", 
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

  fprintf (logFile, "corruptSamples with pSample %f, pChange %f changes %.4f of values\n",
	   pSample, pChange, nChange/(N*(double)M)) ;
  
  pNew->sites = pOld->sites ; pOld->sites = 0 ; 
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(x) ; free (isCorrupt) ;
  return pNew ;
}

PBWT *pbwtCopySamples (PBWT *pOld, int Mnew, double meanLength)
{
  if (!pOld || !pOld->yz) die ("copySample called without an existing pbwt") ;
  PBWT *pNew = pbwtCreate (Mnew, pOld->N) ;
  if (meanLength < 1.0) die ("meanLength %f must be > 1 in pbwtCopySample", meanLength) ;
  int rSwitch = RAND_MAX/meanLength ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;
  uchar *xOld = myalloc (pOld->M, uchar) ;
  int *copy  = myalloc (pNew->M, int) ; /* which old haplotype to copy from */
  int i, j, k ;

  for (j = 0 ; j < pNew->M ; ++j) copy[j] = rand() % pOld->M ; /* initialise copying */

  for (k = 0 ; k < pOld->N ; ++k)
    { for (i = 0 ; i < pOld->M ; ++i) xOld[uOld->a[i]] = uOld->y[i] ;
      for (j = 0 ; j < pNew->M ; ++j) 
	if (rand() < rSwitch) copy[j] = rand() % pOld->M ; /* switch copying */
      for (j = 0 ; j < pNew->M ; ++j) uNew->y[j] = xOld[copy[uNew->a[j]]] ;
      pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uOld) ;
    }  
  pbwtCursorToAFend (uNew, pNew) ;

  fprintf (logFile, "copySamples made %d samples with mean switch length %.1f\n",
	   Mnew, meanLength) ;
  
  pNew->sites = pOld->sites ; pOld->sites = 0 ; 
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->samples = pOld->samples ; pOld->samples = 0 ;
  pbwtDestroy (pOld) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  free(xOld) ; free (copy) ;
  return pNew ;
}

/********** functions to store and retrieve packed dosages ************/

/* these are actually posterior probabilities per haplotype
   to get dosages for a genotype, add the two values, e.g. dg[n] = d[2*n] + d[2*n+1]
   to get genotype likelihoods 
   gl[n][0] = (1-d[2*n]) * (1-d[2*n+1])
   gl[n][1] = d[2*n] + d[2*n+1] - 2*d[2*n]*d[2*n+1]
   gl[n][2] = d[2*n] * d[2*n+1]

   Note that only differences from the {0,1} calls are stored. This is why
   pbwtDosageStore() does not need the cursor but pbwtDosageRetrieve() does.
*/

static inline uchar dosageEncode (double d)
{ if (d > 0.5) d = 1.0 - d ;
  if (!d) return 0 ;
  else return (uchar) (10.0 * (d + 0.0999999)) ; /* value from 0..5 */
}

static inline double dosageDecode (uchar x, uchar y)
{ static double value[16] = { 0.0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.0, 0.0,
			      1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 1.0, 1.0 } ;
  return value[x + (y<<3)] ;
}

static inline uchar *dosageStore (uchar *z, uchar d, int count)
{ if (!d)
    { while (count >= (1 << 15)) { *z++ = 0xff ; count -= 31 << 10 ; }
      if (count >= (1 << 10)) 
	{ *z++ = (7 << 5) | (count >> 10) ; count &= 1023 ; }
      if (count >= (1 << 5)) 
	{ *z++ = (6 << 5) | (count >> 5) ; count &= 31 ; }
      *z++ = count ;
    }
  else
    { while (count >= (1 << 5)) { *z++ = (d << 5) | 31 ; count -= 31 ; }
      *z++ = (d << 5) | count ;
    }
  return z ;
}

void pbwtDosageStore (PBWT *p, double *dosage, int k)
{
  if (!p->dosageOffset) die ("dosageStore called without p->dosageOffset") ;
  long max = arrayMax(p->zDosage) ;
  array(p->dosageOffset,k,long) = max ;
  arrayExtend (p->zDosage, max + p->M) ; /* ensures enough space */
  uchar *z = arrp(p->zDosage, max, uchar) ;
  int i = 0, count = 0 ;
  uchar dLast = 0xff ;
  while (i < p->M)
    { uchar d = dosageEncode (dosage[i]) ;
      if (d != dLast)
	{ if (dLast != 0xff) z = dosageStore (z, dLast, count) ;
	  dLast = d ; count = 0 ;
	}
      ++count ;
      ++i ;
    }
  z = dosageStore (z, dLast, count) ;
  arrayMax(p->zDosage) += z - arrp(p->zDosage, max, uchar) ;
}

double *pbwtDosageRetrieve (PBWT *p, PbwtCursor *u, double *dosage, int k)
{
  if (!dosage) dosage = myalloc (p->M, double) ;
  if (!p->dosageOffset) die ("dosageRetrieve called without p->dosageOffset") ;
  uchar *z = arrp(p->zDosage, arr(p->dosageOffset,k,long), uchar) ;
  int i = 0 ;
  while (i < p->M)
    { uchar x = *z >> 5 ;
      int count = *z & 0x1f ;
      if (x == 6) count <<= 5 ; else if (x == 7) count <<= 10 ;
      while (count--)
	{ /* if (i >= p->M) die ("problem in dosageRetrieve") ; */ 
	  dosage[i] = dosageDecode (x, u->y[i]) ;
	  ++i ;
	}
      ++z ;
    }

  return dosage ;
}

/******************* end of file *******************/

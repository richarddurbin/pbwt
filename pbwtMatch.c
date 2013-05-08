/*  File: pbwtMatch.c
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
 * Description: match functions in pbwt package
 * Exported functions:
 * HISTORY:
 * Last edited: May  5 15:08 2013 (rd)
 * Created: Thu Apr  4 11:55:48 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

/************ finding long (or maximal) matches within the set ************/

static Array matchLengthHist = 0 ;
static uchar **checkHap ;	/* section global for checking only */
static int Ncheck ;		/* section global for checking only */

static void checkMatchMaximal (uchar *x, uchar *y, int start, int end, int N)
{
  int i ;

  if (start && (x[start-1] == y[start-1]))
    die ("match not maximal - can extend backwards\n") ;
  if (end < N && (x[end] == y[end]))
    die ("match not maximal - can extend forwards\n") ;
  for (i = start ; i < end ; ++i)
    if (x[i] != y[i])
      die ("match not a match at %d\n", i) ;
}

static void reportMatch (int ai, int bi, int start, int end)
{
  if (start == end) return ;
  printf ("MATCH\t%d\t%d\t%d\t%d\t%d\n", ai, bi, start, end, end-start) ;

  if (isCheck)			/* check match is a real match and maximal */
    checkMatchMaximal (checkHap[ai], checkHap[bi], start, end, Ncheck) ;
}

static void reportLongMatches1 (Update *u, int T, int k) /* algorithm 3 */
{
  static int *a, *b = 0 ;
  int i, ia, ib, na = 0, nb = 0 ;

  if (!b)			/* initialise */
    { a = myalloc (u->M, int) ;
      b = myalloc (u->M, int) ;
    }
  for (i = 0 ; i < u->M ; ++i)
    { if (u->d[i] > T)
	{ for (ia = 0 ; ia < na ; ++ia)
	    for (ib = 0 ; ib < nb ; ++ib) 
	      reportMatch (u->a[ia], u->a[ib], 0, k) ; /* 0 is wrong! - can't get start */
	  na = 0 ; nb = 0 ;		               /* NB because of this matches won't check */
	}
      if (u->y[i] == 0)
	a[na++] = u->a[i] ;
      else
	b[nb++] = u->a[i] ;
    }
}

static void reportLongMatches2 (Update *u, int T, int k, BOOL isInternal) 
/* alternative giving start - it turns out in tests that this is also faster, so use it */
{
  int i, i0 = 0, ia, ib, na = 0, nb = 0, dmin ;

  for (i = 0 ; i < u->M ; ++i)
    { if (u->d[i] > T)
	{ if (na && nb)		/* then there is something to report */
	    for (ia = i0 ; ia < i ; ++ia)
	      for (ib = ia+1, dmin = 0 ; ib < i ; ++ib)
		{ if (u->d[ib] > dmin) dmin = u->d[ib] ;
		  if (u->y[ib] != u->y[ia])
		    reportMatch (u->a[ia], u->a[ib], dmin, k) ;
		}
	  na = 0 ; nb = 0 ; i0 = i ;
	}
      if (u->y[i] == 0)
	na++ ;
      else
	nb++ ;
    }
}

static void reportMaximalMatches1 (Update *u, int k, BOOL isInternal)
/* algorithm 4 in paper */
{
  int i, j, m, n ;

  for (i = 0 ; i < u->M ; ++i)
    { m = i-1 ; n = i+1 ;
      if (u->d[i] <= u->d[i+1])
	while (u->d[m+1] <= u->d[i])
	  if (u->y[m--] == u->y[i] && isInternal) goto nexti ; /* ERROR in ISMB submission: != should be == */
      if (u->d[i] >= u->d[i+1])
	while (u->d[n] <= u->d[i+1])
	  if (u->y[n++] == u->y[i] && isInternal) goto nexti ; /* ERROR in ISMB submission: != should be == */
      if (isStats)
	++array(matchLengthHist, (u->d[i]<u->d[i+1]) ? k-u->d[i] : k-u->d[i+1], int) ;
      else
	{ for (j = m+1 ; j < i ; ++j) reportMatch (u->a[i], u->a[j], u->d[i], k) ;
	  for (j = i+1 ; j < n ; ++j) reportMatch (u->a[i], u->a[j], u->d[i+1], k) ;
	}
      nexti: ;
    }
}

/* I think there is a good alternative, where I just go down through the list, keeping
   track of some things, but it is not worked out yet, let alone implemented.
*/

void pbwtLongMatches (PBWT *p, int L) /* reporting threshold L - if 0 then maximal */
{
  int k, n = 0 ;
  Update *u = updateCreate (p->M, 0) ;

  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  if (L <= 0) die ("L %d for longWithin must be >0", L) ;

  if (isCheck) { checkHap = pbwtHaplotypes(p) ; Ncheck = p->N ; }

  if (isStats)
    matchLengthHist = arrayReCreate (matchLengthHist, 1000000, int) ;

  for (k = 0 ; k < p->N ; ++k)
    { n += unpack3 (arrp(p->yz,n,uchar), p->M, u->y, 0) ;
      if (L)
	reportLongMatches2 (u, k-L, k, TRUE) ;
      else
	reportMaximalMatches1 (u, k, TRUE) ;
      updateForwardsADU (u, k) ;
    }
  if (L)
    reportLongMatches2 (u, p->N - L, p->N, FALSE) ;
  else
    reportMaximalMatches1 (u, p->N, FALSE) ;

  if (isStats)
    { int i, nTot = 0 ;
      long hTot = 0 ;

      for (i = 0 ; i < arrayMax(matchLengthHist) ; ++i)
	if (arr(matchLengthHist, i, int))
	  { nTot += arr(matchLengthHist, i, int) ;
	    hTot += arr(matchLengthHist, i, int) * i ;
	    printf ("%d\t%d\n", i, arr(matchLengthHist, i, int)) ;
	  }
      fprintf (stderr, "Average %.1f matches per sample\n", nTot/(double)p->M) ;
      fprintf (stderr, "Average length %.1f\n", hTot/(double)nTot) ;
    }

  updateDestroy (u) ;
}

/***************** match new sequences into a reference PBWT ******************/
/***************** this is quite complicated - so multiple implementations ****/

/* First implementation is without using PBWT - works on query and reference haplotypes.
   It should be O(NMQ) in time, O(NM) memory, ignoring query memory.
   NB This version only gives a representative maximal match for each query, start and end
   (the lowest reference index value).
*/

void matchSequencesNaive (PBWT *p, FILE *fp) /* fp is a pbwt file of sequences to match */
{
  PBWT *q = pbwtRead (fp) ;	/* q for "query" of course */
  uchar **query = pbwtHaplotypes (q) ; /* make the query sequences */
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference */
  uchar *x, *y ;		      /* use for current query, y for current reference */
  int i, j, k, kLastMismatch, kk, iBest, N = p->N ;
  int *bestEnd = myalloc (N+1, int) ;
  int *bestSeq = mycalloc (N+1, int) ;
  int totLen = 0, nTot = 0 ;

  if (q->N != p->N) die ("query length in matchSequences %d != PBWT length %d", q->N, p->N) ;

  fprintf (stderr, "Made haplotypes: ") ; timeUpdate () ;

  /* go query by query */
  for (j = 0 ; j < q->M ; ++j)
    { x = query[j] ;
      /* for each query go through targets in turn, updating bestMatch if available */
      memset (bestEnd, 0, N*sizeof(int)) ;
      bestEnd[N] = N+1 ;
      for (i = 0 ; i < p->M ; ++i)
	{ y = reference[i] ;
	  kLastMismatch = N ;	/* NB - never bigger than bestEnd[N] */
	  for (k = p->N ; k-- ; ) /* run backwards so we know end of match */
	    if (x[k] != y[k])
	      { if (kLastMismatch > bestEnd[k+1])
		  for (kk = k+1 ; bestEnd[kk] <= kLastMismatch ; ++kk)
		    { bestEnd[kk] = kLastMismatch ; bestSeq[kk] = i ; }
		kLastMismatch = k ; 
	      }
	  if (kLastMismatch > bestEnd[0])
	    for (kk = 0 ; bestEnd[kk] <= kLastMismatch ; ++kk)
	      { bestEnd[kk] = kLastMismatch ; bestSeq[kk] = i ; }
	}
      /* report best match and update stats */
      int iBest = p->M ;
      for (k = 0 ; k < p->N ; ++k)
	if (bestSeq[k] != iBest)
	  { iBest = bestSeq[k] ;
	    printf ("MATCH query %d to reference %d from %d to %d length %d\n",
		    j, iBest, k, bestEnd[k], bestEnd[k] - k) ;
	    ++nTot ; totLen += bestEnd[k] - k ;
	    if (isCheck)
	      checkMatchMaximal (query[j], reference[iBest], k, bestEnd[k], N) ;
	  }
    }

  fprintf (stderr, "Average number of best matches %.1f, Average length %.1f\n", 
	   nTot/(double)q->M, totLen/(double)nTot) ;

  pbwtDestroy (q) ;
  for (j = 0 ; j < q->M ; ++j) free(query[j]) ; free (query) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  free (bestSeq) ; free (bestEnd) ;
}

/* Next implementation is algorithm 5 from the paper, precalculating indices in memory. 
   It should be O(NQ) time after O(NM) time index calculation. Downside is O(NM) memory,
   13NM bytes for now I think.  This can almost certainly be reduced with some work.
*/

void matchSequencesIndexed (PBWT *p, FILE *fp)
{
  PBWT *q = pbwtRead (fp) ;	/* q for "query" of course */
  uchar **query = pbwtHaplotypes (q) ; /* make the query sequences */
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference */
  uchar *x, *y ;                /* use for current query, and selected reference query */
  Update *up = updateCreate (p->M, 0) ;
  int **a, **d, **u ;		/* stored indexes */
  int e, f, g ;			/* start of match, and pbwt interval as in algorithm 5 */
  int e1, f1, g1 ;		/* next versions of the above, e' etc in algorithm 5 */
  int i, j, k, n = 0, N = p->N, M = p->M ;
  int totLen = 0, nTot = 0 ;

  if (q->N != p->N) die ("query length in matchSequences %d != PBWT length %d", q->N, p->N) ;

  /* build indexes */

  a = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) a[i] = myalloc (p->M, int) ;
  d = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) d[i] = myalloc (p->M+1, int) ;
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  p->c = myalloc (p->N, int) ;
  for (k = 0 ; k < N ; ++k)
    { n += unpack3 (arrp(p->yz,n,uchar), M, up->y, &p->c[k]) ;
      memcpy (a[k], up->a, M*sizeof(int)) ;
      memcpy (d[k], up->d, (M+1)*sizeof(int)) ;
      updateForwardsADU (up, k) ;
      memcpy (u[k], up->u, (M+1)*sizeof(int)) ;
    }
  memcpy (a[k], up->a, M*sizeof(int)) ;
  memcpy (d[k], up->d, (M+1)*sizeof(int)) ;
  updateDestroy (up) ;

  fprintf (stderr, "Made haplotypes and indices: ") ; timeUpdate () ;

  /* match each query in turn */

  for (j = 0 ; j < q->M ; ++j)
    { x = query[j] ;
      e = 0 ; f = 0 ; g = M ;
      for (k = 0 ; k < N ; ++k)
	{               /* use classic FM updates to extend [f,g) interval to next position */
	  f1 = x[k] ? p->c[k] + (f - u[k][f]) : u[k][f] ;
	  g1 = x[k] ? p->c[k] + (g - u[k][g]) : u[k][g] ; 
	  		/* if the interval is non-zero we can just proceed */
	  if (g1 > f1)
	    { f = f1 ; g = g1 ; } /* no change to e */
	  else		/* we have reached a maximum - need to report and update e, f*,g* */
	    { for (i = f ; i < g ; ++i)		/* first report matches */
		{ printf ("MATCH query %d to reference %d from %d to %d length %d\n",
			  j, a[k][i], e, k, k-e) ;
		  if (isCheck) checkMatchMaximal (x, reference[a[k][i]], e, k, N) ;
		}
	      ++nTot ; totLen += k-e ;
	      		/* then update e,f,g */
	      e1 = d[k+1][f1] - 1 ; /* y[f1] and y[f1-1] diverge here, so upper bound for e */
	      if ((x[e1] == 0 && f1 > 0) || f1 == M)
		{ f1 = g1 - 1 ;
		  y = reference[a[k+1][f1]] ; while (x[e1-1] == y[e1-1]) --e1 ;
		  while (d[k+1][f1] <= e1) --f1 ;
		}
	      else if (f1 < M)
		{ g1 = f1 + 1 ;
		  y = reference[a[k+1][f1]] ; while (x[e1-1] == y[e1-1]) --e1 ;
		  while (g1 < M && d[k+1][g1] <= e1) ++g1 ;
		}
	      e = e1 ; f = f1 ; g = g1 ;
	    }
	}
      /* report the maximal matches to the end */
      for (i = f ; i < g ; ++i)
	{ printf ("MATCH query %d to reference %d from %d to %d length %d\n",
		  j, a[k][i], e, k, k-e) ;
	  if (isCheck) checkMatchMaximal (x, reference[a[k][i]], e, k, N) ;
	}
      ++nTot ; totLen += k-e ;
    }

  fprintf (stderr, "Average number of best matches %.1f, Average length %.1f\n", 
	   nTot/(double)q->M, totLen/(double)nTot) ;

  /* cleanup */
  for (j = 0 ; j < q->M ; ++j) free(query[j]) ; free (query) ;
  pbwtDestroy (q) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  for (j = 0 ; j < N ; ++j) free(a[j]) ; free (a) ;
  for (j = 0 ; j < N ; ++j) free(d[j]) ; free (d) ;
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

/* Next is also based on algorithm 5, but applied in parallel to a set of sequences, and
   calculating indices on the fly from PBWT, so low memory for arbitrary reference size.
   It should be O(N(M+Q)) time but only O(N+M) memory.
*/

typedef struct {
  uchar *x ;
  int e, f, g ;			/* start of match, PBWT interval */
#ifdef USE_REVERSE
  int nz, f2, g2 ;		/* pos in reverse pbwt ->zz at e, PBWT interval in this  */
#endif
} MatchInfo ;

void matchSequencesDynamic (PBWT *p, FILE *fp)
/* Rather than keep a[][], d[][], u[][] in memory, which can be very costly, here we update
   them in one sweep through the data. This makes the algorithm O(MN) sadly, but only
   O((M+Q)N) where Q is the number of queries, not O(MNQ) as the naive method.  
       There are actually two versions here, one with isCheck set, in which case the 
   reference haplotypes are calculated, and we can use the algorithm as in Indexed above,
   which requires the check "x[e1-1] == y[e1-1]" when looking for e1.  However this still
   uses MN bytes to store the reference haplotypes.  When isCheck is FALSE, we avoid this
   by extending backwards from k using the PBWT which has some sort of NQlogM cost.
       This could all be made tighter.  First by not explicitly building the u[] array but
   instead updating the query [f,g) in their PBWT sort order.  Second by building a proper
   index on the reference PBWT to allow rapid extension forwards and backwards.  However,
   in order to report which sequence we have hit we need to maintain the a[] array, which
   is O(NM) - this can be done with block memcpy's in the update which might increase
   speed considerably.  Unless we are happy just to find the query maximal strings.
       There are some sections commented out by #ifdef USE_REVERSE.  These show how to keep
   the same interval in the reverse PBWT, as Heng does.  In principle I think we could search
   forwards in this to find e1 when it is closer to the old e than to k, but in practice 
   this is the less common situation, and it seems quite complicated to implement, so is 
   not done.
*/
{
  PBWT *q = pbwtRead (fp) ;	/* q for "query" of course */
  uchar **query = pbwtHaplotypes (q) ; /* make the sequences */
  Array mInfo = arrayCreate (q->M, MatchInfo) ;
  int i, j, k, M = p->M, N = p->N ;
  int e1, f1, g1 ;		/* new values as in algorithm 5 e', f', g' */
  int ny = 0, c ;
  uchar *y ;
  Update *u = updateCreate (M, 0) ;
  MatchInfo *m ;
  uchar **reference ;		/* haplotypes for reference; only made when checking */
  int *oldA = myalloc (M, int) ;
  int totLen = 0, nTot = 0 ;

  if (q->N != p->N) die ("query length in matchSequences %d != PBWT length %d", q->N, p->N) ;

  if (isCheck)
    reference = pbwtHaplotypes (p) ; /* expensive for big reference (same as naive) */

#ifdef USE_REVERSE
  buildReverse (p) ;
#endif

  for (j = 0 ; j < q->M ; ++j)	/* initialise query data structures */
    { m = arrp(mInfo,j,MatchInfo) ;
      m->x = query[j] ;
      m->g = M ;
#ifdef USE_REVERSE
      m->g2 = M ; 
      m->nz = arrayMax(p->zz) ;
#endif
    }
  if (!p->c) p->c = myalloc (N, int) ;

			/* outer loop is now single sweep through pbwt */
  for (k = 0 ; k < N ; ++k)
    { ny += unpack3 (arrp(p->yz,ny,uchar), M, u->y, &p->c[k]) ;
      memcpy (oldA, u->a, M*sizeof(int)) ;
      updateForwardsADU (u, k) ;
      			/* next loop over queries */
      for (j = 0 ; j < q->M ; ++j)
	{ m = arrp(mInfo,j,MatchInfo) ;
			/* use classic FM updates to extend [f,g) interval to next position */
	  f1 = m->x[k] ?  p->c[k] + (m->f - u->u[m->f]) : u->u[m->f] ;
	  g1 = m->x[k] ?  p->c[k] + (m->g - u->u[m->g]) : u->u[m->g] ;
	  		/* if the interval is non-zero we can just proceed */
	  if (g1 > f1)
	    { m->f = f1 ; m->g = g1 ; /* no change to e */
#ifdef USE_REVERSE
	      /* From Heng we can update [f2,g2) because g2-f2 must be the same as g-f */
	      /* because it covers matches to the same substring in reverse, */
	      /* and f2 is unchanged if x[k] == 0, or g2 is unchanged if x[k] == 1 */
	      if (m->x[k]) m->f2 = m->g2 - (g1 - f1) ;
	      else m->g2 = m->f2 + (g1 - f1) ;
#endif
	    }
	  else		/* we have reached a maximum - need to report and update e, f*,g* */
	    { for (i = m->f ; i < m->g ; ++i)		/* first report matches */
		{ printf ("MATCH query %d to reference %d from %d to %d length %d\n",
			  j, oldA[i], m->e, k, k-m->e) ;
		  if (isCheck) checkMatchMaximal (m->x, reference[oldA[i]], m->e, k, N) ;
		}
	      ++nTot ; totLen += k - m->e ;
	      		/* then update e,f,g */
	      e1 = u->d[f1] - 1 ; /* initial upper bound for e1 */
	      if ((m->x[e1] == 0 && f1 > 0) || f1 == M)
		{ f1 = g1 - 1 ;
		  if (isCheck) 
		    { y = reference[u->a[f1]] ; 
		      while (m->x[e1-1] == y[e1-1]) --e1 ;
		    }
		  else		/* have to track back f1 to e1 */
		    { int ny2 = ny, kk = k, f2 = f1 ; uchar z ;
		      do ny2 -= extendPackedBackwards (arrp(p->yz,ny2,uchar), M, &f2, p->c[kk], &z) ;
		      while (m->x[kk--] == z) ;
		      e1 = kk + 2 ;
		    }
		  while (u->d[f1] <= e1) --f1 ;
		}
	      else if (f1 < M)
		{ g1 = f1 + 1 ;
		  if (isCheck) 
		    { y = reference[u->a[f1]] ; 
		      while (m->x[e1-1] == y[e1-1]) --e1 ;
		    }
		  else		/* have to track back f1 to e1 */
		    { int ny2 = ny, kk = k, f2 = f1 ; uchar z = m->x[k] ;
		      do ny2 -= extendPackedBackwards (arrp(p->yz,ny2,uchar), M, &f2, p->c[kk], &z) ;
		      while (m->x[kk--] == z) ;
		      e1 = kk + 2 ;
		    }
		  while (u->d[g1] <= e1) ++g1 ;
		}
	      m->e = e1 ; m->f = f1 ; m->g = g1 ;
	    }
	}
    }

  /* report the maximal matches to the end */
  for (j = 0 ; j < q->M ; ++j)
    { m = arrp(mInfo,j,MatchInfo) ;
      for (i = m->f ; i < m->g ; ++i)
	{ printf ("MATCH query %d to reference %d from %d to %d length %d\n",
		  j, u->a[i], m->e, N, N - m->e) ;
	  if (isCheck) checkMatchMaximal (m->x, reference[u->a[i]], m->e, N, N) ;
	}
      ++nTot ; totLen += N - m->e ;
    }

  fprintf (stderr, "Average number of best matches %.1f, Average length %.1f\n", 
	   nTot/(double)q->M, totLen/(double)nTot) ;

  		/* cleanup */
  for (j = 0 ; j < q->M ; ++j) free(query[j]) ; free (query) ;
  pbwtDestroy (q) ;
  free (oldA) ;
  updateDestroy (u) ;
  if (isCheck) { for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ; }
}

/* Heng Li proposed another approach, based on matching forward and back, but it is intricate
   to program, and this implementation is incomplete.  I am not sure what the complexity is.
*/

void matchSequencesHL (PBWT *p, FILE *fp)
/* ********** INCOMPLETE !! ************ */
{
  PBWT *q = pbwtRead (fp) ;	/* q for "query" of course */
  uchar **query = pbwtHaplotypes (q) ; /* make the sequences */
  uchar *x ;			/* we will use this for the current query */
  int i1 = 0, f1 = 0, g1, nz1 ;		/* [f1,g1) is interval in reverse at i1 */
  int i2, f2, g2, ny2 ;		/* [f2,g2) is interval in forward at i2>i1 */
  int i, j, k, M = p->M, N = p->N ;
  uchar **haps ;		/* haplotypes for target; costly, only made when checking */
  uchar *a = myalloc (M, uchar) ;
  Array info = arrayCreate (64, MatchInfo) ;
  MatchInfo *inf ;

  if (q->N != p->N) die ("query length in matchSequences %d != PBWT length %d", q->N, p->N) ;

  if (isCheck) 
    haps = pbwtHaplotypes (p) ; /* VERY expensive for big targets */
  matchLengthHist = arrayReCreate (matchLengthHist, 100000, int) ;

  if (!p->zz) pbwtBuildReverse (p) ; /* we need the reverse as well as forward BWT */

  for (j = 0 ; j < q->M ; ++j)	/* iterate over queries */
    { x = query[j] ;
      i1 = i2 = 0 ; f1 = f2 = 0 ; g1 = g2 = M ; 
      nz1 = arrayMax(p->zz) ; ny2 = 0 ;
      memcpy (a, p->za, M) ;
      while (i2 < N)
	{ /* first extend i2 forwards until g-f == 0 to find maximal matches */
	  arrayMax(info) = 0 ;
	  while (g2-f2 > 0 && i2 < N)
	    { ny2 += extendMatchForwards (arrp(p->yz,ny2,uchar), M, x[i2], &f2, &g2) ;
	      /* on the way we keep i1,[f1,g1) for whenever g-f drops, in an array */
	      if (g2-f2 < g1-f1)
		{ inf = arrayp(info,arrayMax(info),MatchInfo) ;
		  /* inf->i2 = i2 ; inf->f1 = f1 ; inf->g1 = g1 ; */
		}
	      /* now we can update [f1,g1), knowing g1-f1 must be the same as g2-f2 */
	      /* because it covers matches to the same substring in reverse, */
	      /* and f1 is unchanged if x[i2] == 0, g1 is unchanged if x[i2] == 1 */
	      if (x[i2]) f1 = g1 - (g2 - f2) ;
	      else g1 = f1 + (g2 - f2) ;
	      ++i2 ;
	    }
	  /* now we have to reap maximal matches */
	  /* sweep back from i1 (forward in reverse index), looking at intervals in info  */
	  while (arrayMax(info) && i1 > 0)
	    { inf = arrp(info, arrayMax(info)-1, MatchInfo) ;
	      /*	      for (k = inf->f1 ; k < inf->g1 ; ++k)
		{ ++array(matchLengthHist, inf->i2-i1, int) ;
		  if (isCheck)
		    { printf ("MATCH\t%d\t%d\t%d\t%d\t%d\n", 
			      j, a[k], i1, inf->i2, inf->i2-i1) ;
		      checkMatchMaximal (x, haps[a[k]], i1, inf->i2, N) ;
		    }
		}
	      */
	      if (--arrayMax(info))
		{ 
		}
	    }
	}
    }
  free (a) ;
  for (j = 0 ; j < q->M ; ++j) free(query[j]) ; free (query) ;
  pbwtDestroy (q) ;
  if (isCheck)
    { for (j = 0 ; j < p->M ; ++j) free(haps[j]) ; free (haps) ; }
}

/******************* end of file *******************/

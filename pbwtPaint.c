/*  File: pbwtPaint.c
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
 * Description: tools for chromosome painting as in ChromoPainter, FineStructure etc.
 * Exported functions:
 * HISTORY:
 * Last edited: Apr  3 23:51 2014 (rd)
 * Created: Tue Apr  1 11:34:41 2014 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

#define IDEA3

double **counts = 0 ;

#ifdef IDEA2
static void reportMatch (int i, int j, int start, int end) { ++counts[i][j] ; }
#endif
#ifdef IDEA3
typedef struct { int j ; int start ; int end ; } MatchSegment ;
/* matches are semi-open [start,end) so length is end-start */
static Array *maxMatch = 0 ;	/* arrays of MatchSegment */
static void reportMatch (int i, int j, int start, int end)
{ MatchSegment *ms = arrayp (maxMatch[i], arrayMax(maxMatch[i]), MatchSegment) ;
  ms->j = j ; ms->start = start ; ms->end = end ;
}
#endif

void paintAncestryMatrix (PBWT *p)
{
  int i, j, k ;
  counts = myalloc (p->M, double*) ;
  for (i = 0 ; i < p->M ; ++i)
    counts[i] = mycalloc (p->M, double) ;

#ifdef IDEA1	/* original idea written with Dan Lawson Newton Institute 140402 */
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  int k ;
  for (k = 0 ; k < p->N ; k++)
    { for (j = 1 ; j < p->M ; ++j)
	if (u->y[j] != u->y[j-1])
	  { if (j < p->M || u->d[j] < u->d[j+1]) ++counts[u->a[j]][u->a[j-1]] ;
	    if (j > 1 || u->d[j] < u->d[j-1]) ++counts[u->a[j-1]][u->a[j]] ;
	  }
      pbwtCursorForwardsReadAD (u, k) ;
    }
#endif
#ifdef IDEA2  /* count maximal matches */
  matchMaximalWithin (p, reportMatch) ;
#endif
#ifdef IDEA3  /* weight maximal matches per site as when imputing */
  maxMatch = myalloc (p->M, Array) ;
  for (i = 0 ; i < p->M ; ++i) maxMatch[i] = arrayCreate (1024, MatchSegment) ;
  matchMaximalWithin (p, reportMatch) ;  /* store maximal matches in maxMatch */
  /* now weight per site based on distance from ends */
  for (i = 0 ; i < p->M ; ++i)
    { MatchSegment *m1 = arrp(maxMatch[i],0,MatchSegment), *m ;
      MatchSegment *mStop = arrp(maxMatch[i], arrayMax(maxMatch[i])-1, MatchSegment) ;
      for (k = 1 ; k < p->N ; k++)
	{ double sum = 0 ;
	  while (m1->end <= k && m1 < mStop) ++m1 ;
	  for (m = m1 ; m->start < k && m <= mStop ; ++m) 
	    sum += (k - m->start) * (m->end - k) ;
	  if (sum)
	    for (m = m1 ; m->start < k && m <= mStop ; ++m) 
	      counts[i][m->j] += (k - m->start) * (m->end - k) / sum ;
	}
    }
#endif

  /* report results */
  double *totCounts = mycalloc (p->M, double) ;
  for (i = 0 ; i < p->M ; ++i)
    { for (j = 0 ; j < p->M ; ++j) 
	{ printf (" %8.4g", counts[i][j]) ; 
	  totCounts[i] += counts[i][j] ; 
	}
      putchar ('\n') ;
      if ((i%2) && p->samples) 
	fprintf (stderr, "%s %8.4g %8.4g\n", 
		 sampleName (arr(p->samples,i-1, int)), totCounts[i-1], totCounts[i]) ;
    }

#define HORRIBLE_HACK
#ifdef HORRIBLE_HACK
  int i1, j1 ;
  for (i = 0 ; i < 5 ; ++i)
    { for (j = 0 ; j < 5 ; ++j)
	{ double sum = 0 ;
	  for (i1 = 40*i ; i1 < 40*(i+1) ; ++i1)
	    for (j1 = 40*j ; j1 < 40*(j+1) ; ++j1)
	      sum += counts[i1][j1] ;
	  if (j != i) fprintf (stderr, " %8.4g", sum / (40*40)) ;
	  else fprintf (stderr, " %8.4g", sum / (40*39)) ;
	}
      fputc ('\n', stderr) ;
    }
 #endif

  /* clean up */
  for (i = 0 ; i < p->M ; ++i) free (counts[i]) ; 
  free (counts) ;
  free (totCounts) ;
}

/* end of file */

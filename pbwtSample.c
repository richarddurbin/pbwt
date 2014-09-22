/*  File: pbwtSample.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013
 *-------------------------------------------------------------------
 * Description: functions for samples and populations
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 22 23:07 2014 (rd)
 * * Sep 22 23:07 2014 (rd): move to 64 bit arrays
 * Created: Sat Nov  2 18:42:07 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

/* globals */

static DICT *sampleDict ;
static DICT *populationDict ;
static Array samples ;

/* functions */

void sampleInit (void)
{
  sampleDict = dictCreate (4096) ;
  populationDict = dictCreate (64) ;
  samples = arrayCreate (4096, Sample) ;
  array(samples,0,Sample).nameD = 0 ; /* so that all read samples are non-zero */
}

void sampleDestroy (void)
{
  if (sampleDict) dictDestroy(sampleDict);
  if (populationDict) dictDestroy(populationDict);
  if (samples) arrayDestroy(samples);
}

int sampleAdd (char *name, char *father, char *mother, char *pop)
{
  int k ;
  if (dictAdd (sampleDict, name, &k))
    arrayp(samples, k, Sample)->nameD = k ;
  return k ;
}

Sample *sample (PBWT *p, int i) 
{
  i = arr(p->samples, i, int) ;
  if (i >= arrayMax(samples))
    die ("sample index %d out of range %ld", i, arrayMax(samples)) ;
  return arrp(samples,i,Sample) ;
}

char* sampleName (Sample *s) { return dictName (sampleDict, s->nameD) ; }

char* popName (Sample *s) { return dictName (populationDict, s->popD) ; }

PBWT *pbwtSubSample (PBWT *pOld, Array select)
/* select[i] is the position in old of the i'th position in new */
{
  if (!pOld || !pOld->yz) die ("subSample called without valid pbwt") ;

  PBWT *pNew = pbwtCreate (arrayMax(select), pOld->N) ;
  int i, j, nOld = 0 ;
  uchar *x = myalloc (pNew->M, uchar) ;
  int *ainv = myalloc (pOld->M, int) ;
  PbwtCursor *uOld = pbwtCursorCreate (pOld, TRUE, TRUE) ;
  PbwtCursor *uNew = pbwtCursorCreate (pNew, TRUE, TRUE) ;

  for (i = 0 ; i < pOld->N ; ++i)
    { for (j = 0 ; j < pOld->M ; ++j) ainv[uOld->a[j]] = j ;
      for (j = 0 ; j < pNew->M ; ++j) x[j] = uOld->y[ainv[arr(select,j,int)]] ;
      for (j = 0 ; j < pNew->M ; ++j) uNew->y[j] = x[uNew->a[j]] ;
      pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uOld) ;
    }
  pbwtCursorToAFend (uNew, pNew) ;

  /* need to do this also for missing */

  if (pOld->samples)
    { pNew->samples = arrayCreate (pNew->M, int) ;
      for (j = 0 ; j < pNew->M ; ++j)
	array(pNew->samples,j,int) = arr(pOld->samples,arr(select,j,int),int) ;
    }
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->sites = pOld->sites ; pOld->sites = 0 ;
  pbwtDestroy (pOld) ; /* destroy will free old samples */

  free(x) ; free(ainv) ; pbwtCursorDestroy (uOld) ; pbwtCursorDestroy (uNew) ;
  return pNew ;
}

PBWT *pbwtSubSampleInterval (PBWT *pOld, int start, int Mnew)
{
  if (start < 0 || Mnew <= 0 || start + Mnew > pOld->M)
    die ("bad start %d, Mnew %d in subsample", start, Mnew) ;
  int i ;

  Array select = arrayCreate (Mnew, int) ;
  for (i = 0 ; i < Mnew ; ++i)
    array(select, i, int) = start + i ;

  PBWT *pNew = pbwtSubSample (pOld, select) ;
  arrayDestroy (select) ;
  return pNew ;
}

PBWT *pbwtSelectSamples (PBWT *pOld, FILE *fp)
{
  if (!pOld || !pOld->samples) die ("pbwtSelectSamples called without pre-existing sample names") ;
  int i, j ;

  Array newSamples = pbwtReadSamplesFile (fp) ;

  if (!arrayMax(newSamples)) return pOld ;

  Array oldStart = arrayCreate (arrayMax(samples), int) ; /* start of each sample in old */
  Array oldCount = arrayCreate (arrayMax(samples), int) ; /* how many of each sample in old */
  for (i = 0 ; i < pOld->M ; ++i)
    { if (!array(oldCount, arr(pOld->samples,i,int), int))
	array(oldStart, arr(pOld->samples,i,int), int) = i ;
      ++arr(oldCount, arr(pOld->samples,i,int), int) ;
    }

  Array select = arrayCreate (pOld->M, int) ;
  for (i = 0 ; i < arrayMax(newSamples) ; ++i)
    for (j = 0 ; j < array(oldCount,arr(newSamples,i,int),int) ; ++j)
      array(select,arrayMax(select),int) = arr(oldStart,arr(newSamples,i,int),int)++ ;

  PBWT *pNew = pbwtSubSample (pOld, select) ;
  arrayDestroy (select) ; arrayDestroy (oldCount) ; arrayDestroy (oldStart) ;
  return pNew ;
}

/******************* end of file ****************/

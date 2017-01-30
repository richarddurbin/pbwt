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
static DICT *familyDict ;
static DICT *populationDict ;
static Array samples ;

/* functions */

void sampleInit (void)
{
  sampleDict = dictCreate (4096) ;
  populationDict = dictCreate (64) ;
  familyDict = dictCreate (4096) ;
  samples = arrayCreate (4096, Sample) ;
  array(samples,0,Sample).nameD = 0 ; /* so that all read samples are non-zero */
}

void sampleDestroy (void)
{
  if (sampleDict) dictDestroy(sampleDict);
  if (populationDict) dictDestroy(populationDict);
  if (familyDict) dictDestroy(familyDict);
  if (samples) arrayDestroy(samples);
}

int sampleAdd (char *name)
{ // populate global array
  int k = 0 ;
  if(dictAdd(sampleDict, name, &k)) arrayp(samples, k, Sample)->nameD = k ;
  return k ;
}

Sample *sample (int i) 
{
  if (i >= arrayMax(samples))
    die ("sample index %d out of range %ld", i, arrayMax(samples)) ;
  return arrp(samples,i,Sample) ;
}

int sampleCount (void)
{
    return arrayMax(samples) - 1;
}

int pbwtSamplePloidy(PBWT *p, int i)
{
    if (!p->samples) return 2 ;
    Sample *s = sample (arr(p->samples, i, int)) ;
    if (p->isX) return s->isMale ? 1 : 2 ;
    else if (p->isY) return s->isMale ? 1 : 0 ;
    else return 2 ;
}

char* sampleName (Sample *s) { return dictName (sampleDict, s->nameD) ; }
char* popName (Sample *s, char *name)
{
  if (name) {
    int v ;
    dictAdd(populationDict, name, &v) ;
    s->popD = v ;
  }
  return dictName (populationDict, s->popD) ;
}
char* familyName (Sample *s, char *name)
{
  if (name) {
    int v ;
    dictAdd(familyDict, name, &v) ;
    s->family = v ;
  }
  return dictName (familyDict, s->family) ;
}
Sample* mother (Sample *s, char *name)
{
  if (name) {
    int v = sampleAdd(name) ;
    s->mother = v ;
    s->isFemale = TRUE ;
  }
  return arrp(samples,s->mother,Sample) ;
}
Sample* father (Sample *s, char *name)
{
  if (name) {
    int v = sampleAdd(name) ;
    s->father = v ;
    s->isMale = TRUE ;
  }
  return arrp(samples,s->father,Sample) ;
}

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

  int nMissingSites = 0 ; 
  uchar *missing = mycalloc (pOld->M, uchar) ;
  uchar *xMissing = myalloc (pNew->M+1, uchar) ;
  xMissing[pNew->M] = Y_SENTINEL ;  /* needed for efficient packing */

  for (i = 0 ; i < pOld->N ; ++i)
    { for (j = 0 ; j < pOld->M ; ++j) ainv[uOld->a[j]] = j ;
      for (j = 0 ; j < pNew->M ; ++j) x[j] = uOld->y[ainv[arr(select,j,int)]] ;
      for (j = 0 ; j < pNew->M ; ++j) uNew->y[j] = x[uNew->a[j]] ;
      pbwtCursorWriteForwards (uNew) ;
      pbwtCursorForwardsRead (uOld) ;

      if (pOld->missingOffset)
        { if (!pNew->missingOffset)
            { pNew->zMissing = arrayCreate (10000, uchar) ;
              array(pNew->zMissing, 0, uchar) = 0 ; /* needed so missing[] has offset > 0 */
              pNew->missingOffset = arrayCreate (1024, long) ;
            }
          if (arr(pOld->missingOffset,i,long))
            { 
              unpack3 (arrp(pOld->zMissing, arr(pOld->missingOffset,i,long), uchar), pNew->M, missing, 0) ;
              int nMissing = 0 ;
              for (j = 0 ; j < pNew->M ; ++j) { nMissing += xMissing[j] = missing[arr(select,j,int)] ; }
              if (nMissing) {
                  array(pNew->missingOffset,i,long) = arrayMax(pNew->zMissing) ;
                  pack3arrayAdd (xMissing,pNew->M,pNew->zMissing) ; /* NB original order, not pbwt sort */
                  nMissingSites++ ;
                }
              else
                  array(pNew->missingOffset,pNew->N,long) = 0 ;
            }
          else
              array(pNew->missingOffset,pNew->N,long) = 0 ;
        }
    }
  pbwtCursorToAFend (uNew, pNew) ;
  free(missing) ; free(xMissing) ;

  fprintf (logFile, "%d haplotypes selected from %d, %d missing sites, pbwt size for %d sites is %ld\n", 
           pNew->M, pOld->M, nMissingSites, pNew->N, arrayMax(pNew->yz)) ;

  // TODO: update dosage

  if (pOld->samples)
    { pNew->samples = arrayCreate (pNew->M, int) ;
      for (j = 0 ; j < pNew->M ; ++j)
	array(pNew->samples,j,int) = arr(pOld->samples,arr(select,j,int),int) ;
    }
  pNew->chrom = pOld->chrom ; pOld->chrom = 0 ;
  pNew->sites = pOld->sites ; pOld->sites = 0 ;
  pNew->isX = pOld->isX ; pNew->isY = pOld->isY ;
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

PBWT *pbwtRemoveSamples (PBWT *pOld, FILE *fp)
{
  if (!pOld || !pOld->samples) die ("pbwtRemoveSamples called without pre-existing sample names") ;
  int i, j, k ;

  Array rmSamples = pbwtReadSamplesFile (fp) ;

  if (!arrayMax(rmSamples)) return pOld ;

  Array select = arrayCreate (pOld->M, int) ;
  HASH rmHash = hashCreate (arrayMax(rmSamples)) ;
  for (i = 0 ; i < arrayMax(rmSamples) ; ++i) {
    hashAdd(rmHash,HASH_INT(arr(rmSamples,i,int))) ;
  }

  for (i = 0 ; i < pOld->M ; ++i) {
    if (!hashFind(rmHash,HASH_INT(arr(pOld->samples,i,int))))
      array(select,arrayMax(select),int) = i ;
  }

  hashDestroy(rmHash) ;

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

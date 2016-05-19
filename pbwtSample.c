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
static DICT *keyDict ;
static DICT *valueDict ;
static Array samples ;
static Array metaData ;

/* functions */

void sampleInit (void)
{
  sampleDict = dictCreate (4096) ;
  populationDict = dictCreate (64) ;
  familyDict = dictCreate (4096) ;
  keyDict = dictCreate (4096) ;
  valueDict = dictCreate (4096) ;
  samples = arrayCreate (4096, Sample) ;
  array(samples,0,Sample).nameD = 0 ; /* so that all read samples are non-zero */
  metaData = arrayCreate (4096, MetaData) ;
  array(metaData,0,MetaData).key = 0 ; /* so that all read metadata are non-zero */
}

void sampleDestroy (void)
{
  if (sampleDict) dictDestroy(sampleDict);
  if (populationDict) dictDestroy(populationDict);
  if (familyDict) dictDestroy(familyDict);
  if (keyDict) dictDestroy(keyDict);
  if (valueDict) dictDestroy(valueDict);
  if (samples) arrayDestroy(samples);
  if (metaData) arrayDestroy(metaData);
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

int pbwtSamplePloidy(PBWT *p, int i)
{
    if (!p->samples) return 2 ;
    Sample *s = sample (arr(p->samples, i, int)) ;
    if (p->isX) return s->isMale ? 1 : 2 ;
    else if (p->isY) return s->isMale ? 1 : 0 ;
    else return 2 ;
}

char* sampleName (Sample *s) { return dictName (sampleDict, s->nameD) ; }
Sample *mother (Sample *s) { return arrp(samples,s->mother,Sample) ; }
Sample *father (Sample *s) { return arrp(samples,s->father,Sample) ; }
char* popName (Sample *s) { return dictName (populationDict, s->popD) ; }
char* familyName (Sample *s) { return dictName (familyDict, s->family) ; }

int addMetaData (int sampleID, char *key, char *value, char type)
{
  int k, v ;
  Sample *s = sample (sampleID) ;

  dictAdd(keyDict, key, &k) ;
  int m = arrayMax(metaData) ;
  arrayp(metaData, m, MetaData)->key = k ;
  MetaData *meta = arrp(metaData, m, MetaData) ;

  if (type == 'i') { meta->type = FMF_INT, meta->value.i = strtol(value, NULL, 0); }
  else if (type == 'f') { meta->type = FMF_REAL, meta->value.r = strtod(value, NULL) ; }
  else if (type == 'Z')
  {
    dictAdd(valueDict, value, &v) ;
    meta->type = FMF_STR ; meta->value.s = v ;
    if(!strcasecmp(key, "father")) { dictAdd(sampleDict, value, &v); s->father = v ; }
    if(!strcasecmp(key, "mother")) { dictAdd(sampleDict, value, &v); s->mother = v ; }
    if(!strcasecmp(key, "family")||!strcasecmp(key, "familyID")) { dictAdd(familyDict, value, &v) ; s->family = v ; }
    if(!strcasecmp(key, "pop")||!strcasecmp(key, "population")) { dictAdd(populationDict, value, &v) ; s->popD = v ; }
    if(!strcasecmp(key, "gender")||!strcasecmp(key, "sex"))
    {
      if(!strcasecmp(value, "M")||!strcasecmp(value, "male")) s->isMale = TRUE ;
      if(!strcasecmp(value, "F")||!strcasecmp(value, "female")) s->isFemale = TRUE ;
    }
  }
  else meta->type = FMF_FLAG ;

  if (!s->metaData) s->metaData = arrayCreate (4096, int) ;
  array(s->metaData, arrayMax(s->metaData), int) = m ;

  return m ;
}

MetaData *sampleMetaData (int i) 
{
  if (i >= arrayMax(metaData))
    die ("metaData index %d out of range %ld", i, arrayMax(metaData)) ;
  return arrp(metaData,i,MetaData) ;
}
char *metaDataKey (MetaData *m) { return dictName (keyDict, m->key) ; }
char *metaDataValue (MetaData *m) { return dictName (valueDict, m->value.s) ; }

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

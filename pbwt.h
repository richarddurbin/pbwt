/*  File: pbwt.h
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
 * Description: header file for pbwt package
 * Exported functions:
 * HISTORY:
 * Last edited: May  6 19:41 2013 (rd)
 * Created: Thu Apr  4 11:02:39 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"

/* data types */

typedef unsigned char uchar ;

typedef struct PBWTstruct {
  int N ;			/* number of sites */
  int M ;			/* number of samples */
  char* chrom ;			/* chromosome name */
  Array sites ;			/* array of Site */
  DICT *variationDict ;		/* "xxx|yyy" where variation is from xxx to yyy in VCF */
  Array samples ;		/* array of Sample */
  DICT *sampleNameDict ;	/* as it says */
  Array yz ;			/* compressed PBWT array of uchar */
  Array zz ;			/* compressed reverse PBWT array of uchar */
  int *za ;			/* this is the alphabetical sort order for the strings - made in buildReverse */
  int *c ;			/* c[i] is number of 0s at site i */
} PBWT ;

/* philosophy is to be lazy about PBWT - only fill items for which we have info */
/* NB using a DICT for variation means that identical variations use the same string */

typedef struct SiteStruct {
  int x ;			/* position on chromosome */
  int varD ;			/* index in variationDict */
  int f ;			/* number of 1s, < M */
} Site ;

typedef struct SampleStruct {
  int nameD ;			/* index in sampleNameDict */
} Sample ;

typedef struct {		/* data structure for moving forwards - doesn't know PBWT */
  int M ;
  uchar *y ;
  int *a ;
  int *d ;
  int *u ;
  int *b ;			/* for local operations - no long term meaning */
  int *e ;			/* for local operations - no long term meaning */
} Update ;

/* pbwt_core.c */

extern BOOL isCheck ;		/* when TRUE carry out various checks */
extern BOOL isStats ;		/* when TRUE report stats in various places */

void pbwtInit (void) ;
PBWT *pbwtCreate (int M) ;
void pbwtDestroy (PBWT *p) ;
PBWT *pbwtSubSample (PBWT *pOld, int start, int Mnew) ;
PBWT *pbwtSubSites (PBWT *pOld, double fmin, double frac) ;
PBWT *pbwtSubRange (PBWT *pOld, int start, int end) ;
void pbwtBuildReverse (PBWT *p) ;
uchar **pbwtHaplotypes (PBWT *p) ;

	/* operations to move forwards and backwards in the pbwt using the update structure */

void updateInitialise (Update *u) ;
Update *updateCreate (int M, int *aInit) ;
void updateDestroy (Update *u) ;
void updateForwardsA (Update *u) ; /* algorithm 1 in the manuscript */
void updateBackwardsA (Update *u, int c) ; /* undo algorithm 1 */
void updateForwardsADU (Update *u, int k) ; /* algorithm 2 in the manuscript */

	/* low level operations on packed PBWT, argument yzp in these calls */

#define Y_SENTINEL 2			   /* needed to pack efficiently */
int pack3 (uchar *yp, int M, uchar *yzp) ; /* pack M values from yp into yzp */
int unpack3 (uchar *yzp, int M, uchar *yp, int *n0) ; /* unpack M values from yzp into yp, return number of bytes used from yzp, if (n0) write number of 0s into *n0 */
int packCountReverse (uchar *yzp, int M) ; /* return number of bytes to reverse one position */
int extendMatchForwards (uchar *yzp, int M, uchar x, int *f, int *g) ; /* move hit interval f,g) forwards one position, matching x */
int extendPackedForwards (uchar *yzp, int M, int *f, int c) ; /* move f forwards one position - c is number of 0s */
int extendPackedBackwards (uchar *yzp, int M, int *f, int c, uchar *zp) ; /* move f backwards one position - write value into *zp if zp non-zero */

/* pbwt_io.c */

extern int nCheckPoint ;	/* if set non-zero write pbwt and sites files every n sites when parsing external files */

void pbwtWrite (PBWT *p, FILE *fp) ; /* just writes packed PBWT p->yz */
PBWT *pbwtRead (FILE *fp) ;
void pbwtWriteSites (PBWT *p, FILE *fp) ;
void pbwtReadSites (PBWT *p, FILE *fp) ;
void pbwtReadSamples (PBWT *p, FILE *fp) ;
PBWT *pbwtReadMacs (FILE *fp) ;
PBWT *pbwtReadVcfq (FILE *fp) ;	/* reduced VCF style file made by vcf query */
void pbwtWriteHaplotypes (FILE *fp, PBWT *p) ;

/* pbwt_match.c */

void pbwtLongMatches (PBWT *p, int L) ; /* internal matches longer than L, maximal if L=0 */
void matchSequencesNaive (PBWT *p, FILE *fp) ; /* fp is a pbwt file of sequences to match */
void matchSequencesIndexed (PBWT *p, FILE *fp) ;
void matchSequencesDynamic (PBWT *p, FILE *fp) ;
void matchSequencesHL (PBWT *p, FILE *fp) ; /* incomplete */

/* pbwt_impute.c */

void imputeExplore (PBWT *p, int test) ;
PBWT *phase (PBWT *p, int kMethod, int nSparse) ;
PBWT *pbwtCorruptSites (PBWT *pOld, double pSite, double pChange) ;
PBWT *pbwtCorruptSamples (PBWT *pOld, double pSample, double pChange) ;

/******************* end of file *******************/

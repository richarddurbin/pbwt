/*  File: pbwtIO.c
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
 * Description: read/write functions for pbwt package
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 23 20:41 2014 (rd)
 * Created: Thu Apr  4 11:42:08 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"
#include <ctype.h>

int nCheckPoint = 0 ;	/* if set non-zero write pbwt and sites files every n sites when parsing external files */

BOOL isWriteImputeRef = FALSE ;	/* modifies WriteSites() and WriteHaplotypes() for impute  */

/* basic function to store packed PBWT */

void pbwtWrite (PBWT *p, FILE *fp) /* just writes compressed pbwt in yz */
{
  if (!p || !p->yz) die ("pbwtWrite called without a valid pbwt") ;
  if (!p->aFstart || !p->aFend) die ("pbwtWrite called without start and end indexes") ;

  if (fwrite ("PBW2", 1, 4, fp) != 4) /* version 2, with start and end indexes */
    die ("error writing PBWT in pbwtWrite") ;
  if (fwrite (&p->M, sizeof(int), 1, fp) != 1)
    die ("error writing M in pbwtWrite") ;
  if (fwrite (&p->N, sizeof(int), 1, fp) != 1)
    die ("error writing N in pbwtWrite") ;
  if (fwrite (p->aFstart, sizeof(int), p->M, fp) != p->M)
    die ("error writing aFstart in pbwtWrite") ;
  if (fwrite (p->aFend, sizeof(int), p->M, fp) != p->M)
    die ("error writing aFend in pbwtWrite") ;
  int n = arrayMax(p->yz) ;
  if (fwrite (&n, sizeof(int), 1, fp) != 1)
    die ("error writing n in pbwtWrite") ;
  if (fwrite (arrp(p->yz, 0, uchar), sizeof(uchar), arrayMax(p->yz), fp) != arrayMax(p->yz))
    die ("error writing data in pbwtWrite") ;

  fprintf (stderr, "written %d chars pbwt: M, N are %d, %d\n", arrayMax(p->yz), p->M, p->N) ;
}

void pbwtWriteSites (PBWT *p, FILE *fp)
{
  if (!p || !p->sites) die ("pbwtWriteSites called without sites") ;

  int i ;
  for (i = 0 ; i < p->N ; ++i)
    { Site *s = arrp(p->sites, i, Site) ;
      if (isWriteImputeRef)
	fprintf (fp, "site%d\t%d", i+1, s->x) ;
      else
	fprintf (fp, "%s\t%d", p->chrom ? p->chrom : ".", s->x) ;
      fprintf (fp, "\t%s", dictName (variationDict, s->varD)) ;
      fputc ('\n', fp) ;
    }
  if (ferror (fp)) die ("error writing sites file") ;

  fprintf (stderr, "written %d sites from %d to %d\n", p->N, 
	   arrp(p->sites, 0, Site)->x, arrp(p->sites, p->N-1, Site)->x) ;
}

void pbwtWriteSamples (PBWT *p, FILE *fp)
{
  if (!p || !p->samples) die ("pbwtWriteSamples called without samples") ;

  int i ;
  for (i = 0 ; i < p->M ; i += 2) /* assume diploid for now */
    { Sample *s = sample (p, i) ;
      fprintf (fp, "%s", sampleName(s)) ;
      if (s->popD) fprintf (fp, "\tPOP:%s", popName(s)) ;
      if (s->mother) fprintf (fp, "\tMOTHER:%s", sampleName(s)) ;
      if (s->popD) fprintf (fp, "\tFATHER:%s", popName(s)) ;
      fputc ('\n', fp) ;
    }     
  if (ferror (fp)) die ("error writing samples file") ;

  fprintf (stderr, "written %d samples\n", p->M/2) ;
}

void pbwtWriteMissing (PBWT *p, FILE *fp)
{
  if (!p || !p->zMissing) die ("pbwtWriteMissing called without data on missing") ;

  int n = arrayMax(p->zMissing) ;
  if (fwrite (&n, sizeof(int), 1, fp) != 1)
    die ("error writing n in pbwtWriteMissing") ;
  if (fwrite (arrp(p->zMissing, 0, uchar), sizeof(uchar), n, fp) != n)
    die ("error writing zMissing in pbwtWriteMissing") ;
  if (fwrite (arrp(p->missing, 0, int), sizeof(int), p->N, fp) != p->N)
    die ("error writing missing in pbwtWriteMissing") ;

  fprintf (stderr, "written %d chars compressed missing data\n", n) ;
}

void pbwtWriteReverse (PBWT *p, FILE *fp)
{
  if (!p || !p->zz) die ("pbwtWriteReverse called without reverse pbwt") ;

  Array tz = p->yz ; p->yz = p->zz ;
  int* tstart = p->aFstart ; p->aFstart = p->aRstart ;
  int* tend = p->aFend ; p->aFend = p->aRend ;

  fprintf (stderr, "reverse: ") ; pbwtWrite (p, fp) ;
  
  p->yz = tz ; p->aFstart = tstart ; p->aFend = tend ;
}

#define FOPEN_W(tag)  strcpy (fileNameLeaf, tag) ; if (!(fp = fopen (fileName, "w"))) die ("failed to open %s", fileName)

void pbwtWriteAll (PBWT *p, char *fileNameRoot)
{
  char *fileName = myalloc (strlen (fileNameRoot) + 32, char) ;
  strcpy (fileName, fileNameRoot) ;
  strcat (fileName, ".") ;
  char *fileNameLeaf = fileName + strlen(fileName) ;

  FILE *fp ;
  FOPEN_W("pbwt") ; pbwtWrite (p, fp) ; fclose (fp) ;
  if (p->sites) { FOPEN_W("sites") ; pbwtWriteSites (p, fp) ; fclose (fp) ; }
  if (p->samples) { FOPEN_W("samples") ; pbwtWriteSamples (p, fp) ; fclose (fp) ; }
  if (p->missing) { FOPEN_W("missing") ; pbwtWriteMissing (p, fp) ; fclose (fp) ; }
  if (p->zz) { FOPEN_W("reverse") ; pbwtWriteReverse (p, fp) ; fclose (fp) ; }

  free (fileName) ;
}

void pbwtCheckPoint (PBWT *p)
{
  static BOOL isA = TRUE ;
  char fileNameRoot[20] ;

  sprintf (fileNameRoot, "check_%c", isA ? 'A' : 'B') ;
  pbwtWriteAll (p, fileNameRoot) ;

  isA = !isA ;
}

/*******************************/

PBWT *pbwtRead (FILE *fp) 
{
  int m, n ;
  PBWT *p ;
  static char tag[5] = "test" ;
  int version ;

  if (fread (tag, 1, 4, fp) != 4) die ("failed to read 4 char tag - is file readable?") ;
  if (!strcmp (tag, "PBW2")) version = 2 ; /* current version */
  else if (!strcmp (tag, "PBWT")) version = 1 ; /* without start, end indexes */
  else if (!strcmp (tag, "GBWT")) version = 0 ; /* earliest version */
  else die ("failed to recognise file type %s in pbwtRead - was it written by pbwt?", tag) ;

  if (fread (&m, sizeof(int), 1, fp) != 1) die ("error reading m in pbwtRead") ;
  if (fread (&n, sizeof(int), 1, fp) != 1) die ("error reading n in pbwtRead") ;
  p = pbwtCreate (m, n) ;
  if (version > 1)		/* read aFstart and aFend */
    { p->aFstart = myalloc (m, int) ;
      if (fread (p->aFstart, sizeof(int), m, fp) != m) die ("error reading aFstart in pbwtRead") ;
      p->aFend = myalloc (m, int) ;
      if (fread (p->aFend, sizeof(int), m, fp) != m) die ("error reading aFend in pbwtRead") ;
    }
  else				/* set aFstart to 0..M-1, leave aFend empty */
    { p->aFstart = myalloc (m, int) ;
      int i ; for (i = 0 ; i < m ; ++i) p->aFstart[i] = i ;
    }
  if (fread (&n, sizeof(int), 1, fp) != 1)
    die ("error reading pbwt file") ;
  p->yz = arrayCreate (n, uchar) ;
  array(p->yz, n-1, uchar) = 0 ; /* sets arrayMax */
  if (fread (arrp(p->yz, 0, uchar), sizeof(uchar), n, fp) != n)
    die ("error reading data in pbwt file") ;

  fprintf (stderr, "read pbwt %s file with %d bytes: M, N are %d, %d\n", tag, n, p->M, p->N) ;
  return p ;
}

static BOOL readMatchChrom (char **pChrom, FILE *fp)
{
  char *newChrom = fgetword (fp) ;

  if (strcmp (newChrom, "."))
    if (!*pChrom) 
      *pChrom = strdup (newChrom) ;
    else if (strcmp (newChrom, *pChrom))
      return FALSE ;
  return TRUE ;
}

Array pbwtReadSitesFile (FILE *fp, char **chrom)
{
  char c ;
  Site *s ;
  int line = 1 ;
  Array varTextArray = arrayCreate (256, char) ;
  Array sites = arrayCreate (4096, Site) ;

  while (!feof(fp))
    if (readMatchChrom (chrom, fp))	/* if p->chrom then match, else set if not "." */
      { if (feof(fp)) break ;
	s = arrayp(sites, arrayMax(sites), Site) ;
	s->x = 0 ; while (isdigit(c = fgetc(fp))) s->x = s->x*10 + (c-'0') ;
	if (!feof(fp) && c != '\n')
	  { if (!isspace (c)) die ("bad position line %d in sites file", line) ;
	    while (isspace(c = fgetc(fp)) && c != '\n') ;
	    if (c == '\n') die ("bad end of line at line %d in sites file", line) ;
	    int i = 0 ; array(varTextArray, i++, char) = c ;
	    while ((c = fgetc(fp)) && c != '\n') 
	      array(varTextArray, i++, char) = c ;
	    array(varTextArray, i, char) = 0 ;
	    dictAdd (variationDict, arrayp(varTextArray,0,char), &s->varD) ;
	    while (c != '\n' && !feof(fp)) c = fgetc(fp) ;
	  }
	++line ;
      }
    else if (!feof(fp))
      die ("failed to match chromosome in sites file: line %d", line) ;

  if (ferror (fp)) die ("error reading sites file") ;
  
  fprintf (stderr, "read %d sites on chromosome %s from file\n", arrayMax(sites), *chrom) ;

  arrayDestroy (varTextArray) ;
  return sites ;
}

void pbwtReadSites (PBWT *p, FILE *fp)
{
  if (!p) die ("pbwtReadSites called without a valid pbwt") ;

  p->sites = pbwtReadSitesFile (fp, &p->chrom) ;
  if (arrayMax(p->sites) != p->N)
    die ("sites file contains %d sites not %d as in pbwt", arrayMax(p->sites), p->N) ;
}

Array pbwtReadSamplesFile (FILE *fp) /* for now assume all samples diploid */
/* should add code to this to read father and mother and population
   propose to use IMPUTE2 format for this */
{
  char *name, c ;
  int line = 0 ;
  Array nameArray = arrayCreate (1024, char) ;
  Array samples = arrayCreate (1024, int) ;

  while (!feof(fp))
    { int n = 0 ;
      while ((c = getc(fp)) && !isspace(c) && !feof (fp)) array(nameArray, n++, char) = c ;
      if (feof(fp)) break ;
      if (!n) die ("no name line %d in samples file", arrayMax(samples)+1) ;
      array(nameArray, n++, char) = 0 ;
      /* next bit of code deals with header lines for IMPUTE2 samples file, so can read that */
      if (!strcmp (arrp(nameArray,0,char), "ID_1") && !arrayMax(samples))
	{ while ((c = getc(fp)) && c != '\n' && !feof(fp)) ; /* remove header line */
	  while ((c = getc(fp)) && c != '\n' && !feof(fp)) ; /* and next line of zeroes?? */
	  continue ;
	}
      array(samples,arrayMax(samples),int) = sampleAdd (arrp(nameArray,0,char), 0, 0, 0) ;
      /* now remove the rest of the line, for now */
      while (c != '\n' && !feof(fp)) c = getc(fp) ;
    }
  arrayDestroy (nameArray) ;

  fprintf (stderr, "read %d sample names\n", arrayMax(samples)) ;

  return samples ;
}

void pbwtReadSamples (PBWT *p, FILE *fp)
{
  Array samples = pbwtReadSamplesFile (fp) ;
  if (arrayMax(samples) != p->M/2) 
    die ("wrong number of diploid samples: %d needed", p->M/2) ;
  p->samples = arrayReCreate(p->samples, p->M, int) ;
  int i ; 
  for (i = 0 ; i < arrayMax(samples) ; ++i)
    { array(p->samples, 2*i, int) = arr(samples, i, int) ;
      array(p->samples, 2*i+1, int) = arr(samples, i, int) ;
    }
  arrayDestroy (samples) ;
}

void pbwtReadMissing (PBWT *p, FILE *fp)
{
  if (!p) die ("pbwtReadSamples called without a valid pbwt") ;

  int n ; if (fread (&n, sizeof(int), 1, fp) != 1)
    die ("error reading n in pbwtReadMissing") ;
  p->zMissing = arrayReCreate (p->zMissing, n, uchar) ;
  if (fread (arrp(p->zMissing, 0, uchar), sizeof(uchar), n, fp) != n)
    die ("error reading zMissing in pbwtReadMissing") ;
  array(p->zMissing, n-1, uchar) ; /* forces setting max correctly - not sure if needed */
  p->missing = arrayReCreate (p->missing, p->N, int) ;
  if (fread (arrp(p->missing, 0, int), sizeof(int), p->N, fp) != p->N)
    die ("error reading missing in pbwtReadMissing") ;

  fprintf (stderr, "read %d chars compressed missing data\n", n) ;
}

void pbwtReadReverse (PBWT *p, FILE *fp)
{
  if (!p) die ("pbwtReadReverse called without a valid pbwt") ;

  PBWT *q = pbwtRead (fp) ;
  if (q->M != p->M || q->N != p->N)
    die ("M %d or N %d in reverse don't match %, %d in forward", q->M, q->N, p->M, p->N) ;
  p->zz = q->yz ; q->yz = 0 ;
  p->aRstart = q->aFstart ; q->aFstart = 0 ;
  p->aRend = q->aFend ; q->aFend = 0 ;
  pbwtDestroy (q) ;
 }

#define FOPEN_R(tag)  strcpy (fileNameLeaf, tag), (fp = fopen (fileName, "r"))

PBWT *pbwtReadAll (char *fileNameRoot)
{
  PBWT *p ;

  char *fileName = myalloc (strlen (fileNameRoot) + 32, char) ;
  strcpy (fileName, fileNameRoot) ;
  strcat (fileName, ".") ;
  char *fileNameLeaf = fileName + strlen(fileName) ;

  FILE *fp ;
  if (FOPEN_R ("pbwt")) { p = pbwtRead (fp) ; fclose (fp) ; } 
  else die ("failed to open %s", fileName) ;
  if (FOPEN_R("sites")) { pbwtReadSites (p, fp) ; fclose (fp) ; }
  if (FOPEN_R("samples")) { pbwtReadSamples (p, fp) ; fclose (fp) ; }
  if (FOPEN_R("missing")) { pbwtReadMissing (p, fp) ; fclose (fp) ; }
  if (FOPEN_R("reverse")) { pbwtReadReverse (p, fp) ; fclose (fp) ; }

  free (fileName) ;
  return p ;
}

/****************** MaCS file parsing *****************/

static void parseMacsHeader (FILE *fp, int *M, double *L)	/* parse MACS file header */
{
  if (strcmp (fgetword (fp), "COMMAND:")) die ("MaCS COMMAND line not found") ;
  fgetword (fp) ; 		/* the command */
  *M = atoi(fgetword(fp)) ; if (!*M) die ("failed to get M") ;
  *L = atof(fgetword(fp)) ; if (!*L) die ("failed to get L") ;
  while (fgetc(fp) != '\n') ;	/* ignore rest of line */
  
  if (strcmp (fgetword (fp), "SEED:")) die ("SEED line not found") ;
  while (fgetc(fp) != '\n') ;	/* ignore rest of line */
}

static BOOL parseMacsSite (FILE *fp, Site *s, int M, double L, uchar *yp) /* parse MACS site line */
{
  static uchar conv[256] ;
  static BOOL isFirst = TRUE ;
  int number ;

  if (isFirst) { conv['0'] = 0 ; conv['1'] = 1 ; isFirst = FALSE ; }
  if (feof (fp)) return FALSE ;

  if (strcmp (fgetword (fp), "SITE:")) return FALSE;

  number = atoi(fgetword(fp)) ;	/* this is the site number */
  s->x = (int) (L * atof(fgetword(fp))) ;
  atof(fgetword(fp)) ;		/* ignore the time */
  while (M--)
    *yp++ = conv[getc(fp)] ;
  if (feof (fp)) return FALSE ;
  if (getc(fp) != '\n') die ("end of line error for MaCS SITE %d", number) ;

  return TRUE ;
}

PBWT *pbwtReadMacs (FILE *fp)
{
  PBWT *p ;
  int M ;
  double L ;
  int i,j, nxPack, nyPack, nUnpack, n0 ;
  int nxTot = 0, nyTot = 0 ;
  uchar *x ;
  int *a ;
  PbwtCursor *u ;

  parseMacsHeader (fp, &M, &L) ;
  p = pbwtCreate (M, 0) ;
  p->sites = arrayCreate(4096, Site) ;
  u = pbwtCursorCreate (p, TRUE, TRUE) ;
  x = myalloc (M+1, uchar) ; x[M] = Y_SENTINEL ;	/* sentinel required for packing */

  while (parseMacsSite (fp, arrayp(p->sites,p->N,Site), M, L, x))
    { for (j = 0 ; j < M ; ++j) u->y[j] = x[u->a[j]] ; /* next character in sort order: BWT */
      pbwtCursorWriteForwards (u) ;
      p->N++ ;
      if (nCheckPoint && !(p->N % nCheckPoint))
	pbwtCheckPoint (p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (stderr, "read MaCS file: M, N are\t%d\t%d\n", M, p->N) ;
  if (isStats)
    fprintf (stderr, "                xtot, ytot are\t%d\t%d\n", nxTot, nyTot) ;

  free(x) ; pbwtCursorDestroy (u) ;

  return p ;
}

/************* read vcfq files, made with vcf query ... ***************/

static char* getVariation (FILE *fp)
{
  static Array a = 0 ;
  char c ;
  int i = 0 ;

  if (!a) a = arrayCreate (64, char) ;
  while (!isspace(c = getc(fp)) && !feof(fp)) array(a,i++,char) = c ;
  array(a,i++,char) = c ;
  while (!isspace(c = getc(fp)) && !feof(fp)) array(a,i++,char) = c ;
  array(a,i++,char) = 0 ;

  return arrp(a,0,char) ;
}

static BOOL parseVcfqLine (PBWT **pp, FILE *fp, Array x)
{
  char c, *chrom, *var ;
  int pos, m, len ;
  Site *s ;
  PBWT *p = *pp ;

  if (!p)			/* first call on a file - don't know M yet */
    chrom = strdup (fgetword (fp)) ;
  else
    while (!feof (fp) && !readMatchChrom (&p->chrom, fp))
      while (!feof (fp) && (getc (fp) != '\n')) ; /* ignore rest of line */
  if (feof (fp)) return 0 ;

  pos = atoi (fgetword (fp)) ;
  var = getVariation (fp) ;

  m = 0 ;
  while ((c = getc (fp)))
    switch (c)
      { case EOF: return 0 ;
      case '\n': goto endOfLine ; /* -1 is EOF */
      case '0': array(x,m++,uchar) = 0 ; break ;
      case '1': array(x,m++,uchar) = 1 ; break ;
      case '|': case '/': case '\\': case '\t': break ; /* could check ploidy here */
      default: die ("unexpected character %d in vcfq file genotype section") ;
      }
 endOfLine:
  if (feof (fp)) return FALSE ;
  if (p && m != p->M) die ("length mismatch reading vcfq line") ;

  if (!*pp)
    { p = *pp = pbwtCreate (m, 0) ;
      if (strcmp (chrom, ".")) p->chrom = chrom ;
      p->sites = arrayCreate(4096, Site) ;
      array(x,p->M,uchar) = Y_SENTINEL ; /* sentinel required for packing */
    }
  s = arrayp(p->sites, arrayMax(p->sites), Site) ;
  s->x = pos ;
  dictAdd (variationDict, var, &s->varD) ;

  ++p->N ;
  return TRUE ;
}

typedef BOOL (*ParseLineFunc)(PBWT** pp, FILE *f, Array a) ;

static PBWT *pbwtReadLineFile (FILE *fp, char* type, ParseLineFunc parseLine)
{
  PBWT *p = 0 ;
  int j ;
  uchar *x ;		/* original, sorted, compressed */
  int *a ;
  Array xArray = arrayCreate (10000, uchar) ;
  PbwtCursor *u ;

  while ((*parseLine) (&p, fp, xArray)) /* create p first time round */
    { if (!p->yz)		/* first line; p was just made! */
	{ p->yz = arrayCreate(4096*32, uchar) ;
	  u = pbwtCursorCreate (p, TRUE, TRUE) ;
	}
      x = arrp(xArray,0,uchar) ;
      for (j = 0 ; j < p->M ; ++j) u->y[j] = x[u->a[j]] ;
      pbwtCursorWriteForwards (u) ;
      if (nCheckPoint && !(p->N % nCheckPoint))	pbwtCheckPoint (p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (stderr, "read %s file", type) ;
  if (p->chrom) fprintf (stderr, " for chromosome %s", p->chrom) ;
  fprintf (stderr, ": M, N are\t%d\t%d; yz length is %d\n", p->M, p->N, arrayMax(p->yz)) ;

  arrayDestroy(xArray) ; pbwtCursorDestroy (u) ;

  return p ;
}

PBWT *pbwtReadVcfq (FILE *fp) { return pbwtReadLineFile (fp, "vcfq", parseVcfqLine) ; }

/*************** read impute2 .gen format - contains sites not samples **********/

static long int nGenMissing ;

static BOOL parseGenLine (PBWT **pp, FILE *fp, Array x) /* based on parseVcfqLine() */
{
  PBWT *p = *pp ;

  fgetword(fp) ; fgetword(fp) ;	/* ignore first two name fields */

  int pos = atoi (fgetword (fp)) ;
  char *var = getVariation (fp) ; /* but need to change ' ' separator to '\t' */
  if (feof (fp)) return FALSE ;
  char *cp = var ; while (*cp && *cp != ' ') ++cp ; if (*cp == ' ') *cp = '\t' ; else die ("missing separator in line %d, var is %d", p?p->N:0, var) ;

  int m = 0, nscan ;
  while (!feof(fp))
    { char c = getc(fp) ; if (c == '\n') break ; else if (!isspace(c)) ungetc (c, fp) ;
      float f0, f1, f2 ;
      if ((nscan = fscanf (fp, "%f %f %f", &f0, &f1, &f2) != 3)) 
	die ("bad line %d, %d floats found, m %d, pos %d, var %s", p ? p->N : 0, nscan, m, pos, var) ;
      if (f0 + f1 + f2 == 0)	/* missing genotype */
	{ f0 = 1 ; ++nGenMissing ; }
      if (f0 + f1 + f2 < 0.98) 
	die ("inconsistent genotype in gen file: %f %f %f at %d line %d", 
	     f0, f1, f2, m, p ? p->N : 0) ;
      if (f0 > f1 && f0 > f2) { array(x,m++,uchar) = 0 ; array(x,m++,uchar) = 0 ; }
      else if (f1 > f2) { array(x,m++,uchar) = 0 ; array(x,m++,uchar) = 1 ; }
      else /* f2 is largest */ { array(x,m++,uchar) = 1 ; array(x,m++,uchar) = 1 ; }
    }

  if (feof (fp)) return FALSE ;
  if (p && m != p->M) die ("length mismatch reading vcfq line") ;

  if (!*pp)
    { p = *pp = pbwtCreate (m, 0) ;
      p->sites = arrayCreate(4096, Site) ;
      array(x,p->M,uchar) = Y_SENTINEL ; /* sentinel required for packing */
    }
  Site *s = arrayp(p->sites, arrayMax(p->sites), Site) ;
  s->x = pos ;
  dictAdd (variationDict, var, &s->varD) ;

  ++p->N ;
  return TRUE ;
}

static BOOL parseHapLine (PBWT **pp, FILE *fp, Array x) /* same as parseGenLine - slightly simpler */
{
  PBWT *p = *pp ;
  
  fgetword(fp) ; fgetword(fp) ;	/* ignore first two name fields */
  
  int pos = atoi (fgetword (fp)) ;
  char *var = getVariation (fp) ; /* but need to change ' ' separator to '\t' */
  if (feof (fp)) return FALSE ;
  char *cp = var ; while (*cp && *cp != ' ') ++cp ; if (*cp == ' ') *cp = '\t' ; else die ("missing separator in line %d, var is %d", p?p->N:0, var) ;
  
  int m = 0, nscan ;
  while (!feof(fp))
  { char c = getc(fp) ; if (c == '\n') break ; else if (!isspace(c)) ungetc (c, fp) ;
    float f0, f1 ;
    if ((nscan = fscanf (fp, "%f %f", &f0, &f1) != 2))
      { warn ("bad line %d, %d floats found, m %d, pos %d, var %s - aborting", p ? p->N : 0, nscan, m, pos, var) ;
	return FALSE ;
      }
    array(x,m++,uchar) = f0 ; array(x,m++,uchar) = f1 ; /* haps are phased, put straight in as is */
  }
  
  if (feof (fp)) return FALSE ;
  if (p && m != p->M) die ("length mismatch reading haps line") ;
  
  if (!*pp)
  { p = *pp = pbwtCreate (m, 0) ;
    p->sites = arrayCreate(4096, Site) ;
    array(x,p->M,uchar) = Y_SENTINEL ; /* sentinel required for packing */
  }
  Site *s = arrayp(p->sites, arrayMax(p->sites), Site) ;
  s->x = pos ;
  dictAdd (variationDict, var, &s->varD) ;
  
  ++p->N ;
  return TRUE ;
}


PBWT *pbwtReadGen (FILE *fp, char *chrom) 
{ 
  nGenMissing = 0 ;
  PBWT *p = pbwtReadLineFile (fp, "gen", parseGenLine) ;
  p->chrom = strdup (chrom) ;
  if (nGenMissing) fprintf (stderr, "%ld missing genotypes set to 00\n", nGenMissing) ;
  return p ;
}

PBWT *pbwtReadHap (FILE *fp, char *chrom)
{
  PBWT *p = pbwtReadLineFile (fp, "hap", parseHapLine) ;
  p->chrom = strdup (chrom) ;
  return p ;
}

PBWT *pbwtReadPhase (FILE *fp) /* Li and Stephens PHASE format */
{
  fgetword (fp) ; if (getc(fp) != '\n') die ("bad first line in phase file") ;
  int m = atoi (fgetword (fp)) ;  if (getc(fp) != '\n') die ("bad second line in phase file") ;
  int n = atoi (fgetword (fp)) ;  if (getc(fp) != '\n') die ("bad third line in phase file") ;
  PBWT *p = pbwtCreate (2*m, n) ;
  p->chrom = strdup (fgetword(fp)) ; /* example 4th line is P followed by site positions */
  p->sites = arrayCreate (4096, Site) ;
  int i ; for (i = 0 ; i < p->N ; ++i) arrayp(p->sites,i,Site)->x = atoi (fgetword(fp)) ;
  if (getc(fp) != '\n') die ("bad 4th line in phase file") ;
  char var[2] ; var[1] = 0 ;
  for (i = 0 ; i < p->N ; ++i) 
    { *var = getc(fp) ; dictAdd (variationDict, var, &(arrayp(p->sites,i,Site)->varD)) ; }
  if (getc(fp) != '\n') die ("bad 5th line in phase file") ;
  uchar **data = myalloc (p->N, uchar*) ;
  for (i = 0 ; i < p->N ; ++i) data[i] = myalloc (p->M, uchar) ;
  int j ; for (j = 0 ; j < p->M ; ++j)
    { for (i = 0 ; i < p->N ; ++i) data[i][j] = getc(fp) - '0' ;
      if (getc(fp) != '\n') die ("bad %dth line in phase file", 6+j) ;
    }
  p->yz = arrayCreate(4096*32, uchar) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  for (i = 0 ; i < p->N ; ++i)
    { for (j = 0 ; j < p->M ; ++j) u->y[j] = data[i][u->a[j]] ;
      pbwtCursorWriteForwards (u) ;
      if (nCheckPoint && !((i+1) % nCheckPoint)) pbwtCheckPoint (p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (stderr, "read phase file") ;
  if (p->chrom) fprintf (stderr, " for chromosome %s", p->chrom) ;
  fprintf (stderr, ": M, N are\t%d\t%d; yz length is %d\n", p->M, p->N, arrayMax(p->yz)) ;

  for (i = 0 ; i < p->N ; ++i) free(data[i]) ;
  free (data) ; pbwtCursorDestroy (u) ;
  return p ;
}

/*************** write haplotypes ******************/

void pbwtWriteHaplotypes (FILE *fp, PBWT *p)
{
  int i, j, n = 0, M = p->M ;
  uchar *hap = myalloc (M, uchar) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;

  for (i = 0 ; i < p->N ; ++i)
    { for (j = 0 ; j < M ; ++j) hap[u->a[j]] = u->y[j] ;
      for (j = 0 ; j < M ; ++j) 
	{ if (isWriteImputeRef && j > 0) putc (' ', fp) ;
	  putc (hap[j]?'1':'0', fp) ;
	}
      putc ('\n', fp) ; fflush (fp) ;
      pbwtCursorForwardsRead (u) ;
    }
  free (hap) ; pbwtCursorDestroy (u) ;

  fprintf (stderr, "written haplotype file: %d rows of %d\n", p->N, M) ;
}

/*************** write IMPUTE files ********************/

void pbwtWriteImputeRef (PBWT *p, char *fileNameRoot)
{
  char *fileName = myalloc (strlen (fileNameRoot) + 32, char) ;
  strcpy (fileName, fileNameRoot) ;
  strcat (fileName, ".") ;
  char *fileNameLeaf = fileName + strlen(fileName) ;
  FILE *fp ;

  isWriteImputeRef = TRUE ;

  FOPEN_W("imputeHaps") ; pbwtWriteHaplotypes (fp, p) ; fclose (fp) ;

  FOPEN_W("imputeLegend") ; 
  fprintf (fp, "rsID\tposition\ta0\ta1\n") ; /* header line */
  pbwtWriteSites (p, fp) ; fclose (fp) ;

  isWriteImputeRef = FALSE ;

  free (fileName) ;
}

void pbwtWriteImputeHapsG (PBWT *p, FILE *fp)
{
  if (!p || !p->sites) die ("pbwtWriteImputeHaps called without sites") ;

  int i, j ;
  uchar *hap = myalloc (p->M, uchar) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;

  for (i = 0 ; i < p->N ; ++i)
    { Site *s = arrp(p->sites, i, Site) ;
      fprintf (fp, "site%d\tsite%d\t%d", i+1, i+1, s->x) ;
      fprintf (fp, "\t%s", dictName (variationDict, s->varD)) ;
      
      for (j = 0 ; j < p->M ; ++j) hap[u->a[j]] = u->y[j] ;
      for (j = 0 ; j < p->M ; ++j) putc (' ', fp), putc (hap[j]?'1':'0', fp) ;
      putc ('\n', fp) ; fflush (fp) ;
      pbwtCursorForwardsRead (u) ;
    }

  free (hap) ; pbwtCursorDestroy (u) ;
}

void pbwtWriteGen (PBWT *p, FILE *fp)
{
  if (!p || !p->sites) die ("pbwtWriteImputeHaps called without sites") ;

  int i, j ;
  uchar *hap = myalloc (p->M, uchar) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;

  for (i = 0 ; i < p->N ; ++i)
    { Site *s = arrp(p->sites, i, Site) ;
      fprintf (fp, "site%d\tsite%d\t%d", i+1, i+1, s->x) ;
      fprintf (fp, "\t%s", dictName (variationDict, s->varD)) ;
      
      for (j = 0 ; j < p->M ; ++j) hap[u->a[j]] = u->y[j] ;
      for (j = 0 ; j < p->M ; j+=2) 
	if (hap[j] + hap[j+1] == 0) fprintf (fp, " 1 0 0") ;
	else if (hap[j] + hap[j+1] == 1) fprintf (fp, " 0 1 0") ;
	else fprintf (fp, " 0 0 1") ;
      putc ('\n', fp) ; fflush (fp) ;
      pbwtCursorForwardsRead (u) ;
    }

  free (hap) ; pbwtCursorDestroy (u) ;
}

/******************* end of file *******************/

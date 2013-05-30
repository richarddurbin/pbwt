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
 * Last edited: May  6 23:56 2013 (rd)
 * Created: Thu Apr  4 11:42:08 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

int nCheckPoint = 0 ;	/* if set non-zero write pbwt and sites files every n sites when parsing external files */

/* basic function to store packed PBWT */

void pbwtWrite (PBWT *p, FILE *fp) /* just writes compressed pbwt in yz */
{
  int n ;

  if (!p || !p->yz) die ("pbwtWrite called without a valid pbwt") ;

  n = arrayMax(p->yz) ;
  if (fwrite ("PBWT", 1, 4, fp) != 4)
    die ("error writing PBWT in pbwtWrite") ;
  if (fwrite (&p->M, sizeof(int), 1, fp) != 1)
    die ("error writing M in pbwtWrite") ;
  if (fwrite (&p->N, sizeof(int), 1, fp) != 1)
    die ("error writing N in pbwtWrite") ;
  if (fwrite (&n, sizeof(int), 1, fp) != 1)
    die ("error writing n in pbwtWrite") ;
  if (fwrite (arrp(p->yz, 0, uchar), sizeof(uchar), arrayMax(p->yz), fp) != arrayMax(p->yz))
    die ("error writing data in pbwtWrite") ;

  fprintf (stderr, "written %d chars pbwt: M, N are %d, %d\n", arrayMax(p->yz), p->M, p->N) ;
}

void pbwtWriteSites (PBWT *p, FILE *fp)
{
  int i ;
  Site *s ;

  if (!p || !p->sites) die ("pbwtWriteSites called without a valid pbwt") ;

  for (i = 0 ; i < p->N ; ++i)
    { s = arrp(p->sites, i, Site) ;
      fprintf (fp, "%s\t%d", p->chrom ? p->chrom : ".", s->x) ;
      if (p->variationDict) fprintf (fp, "\t%s", dictName (p->variationDict, s->varD)) ;
      fputc ('\n', fp) ;
    }
  if (ferror (fp)) die ("error writing sites file") ;

  fprintf (stderr, "written %d sites from %d to %d\n", p->N, 
	   arrp(p->sites, 0, Site)->x, arrp(p->sites, p->N-1, Site)->x) ;
}

static void pbwtCheckPoint (PBWT *p)
{
  static BOOL isA = TRUE ;
  char fileName[20] ;
  FILE *fp ;

  sprintf (fileName, "check_%c.pbwt", isA ? 'A' : 'B') ;
  if (!(fp = fopen (fileName, "w"))) die ("failed to open checkpoint file %s", fileName) ;
  pbwtWrite (p, fp) ;
  fclose (fp) ;

  sprintf (fileName, "check_%c.sites", isA ? 'A' : 'B') ;
  if (!(fp = fopen (fileName, "w"))) die ("failed to open checkpoint file %s", fileName) ;
  pbwtWriteSites (p, fp) ;
  fclose (fp) ;

  isA = !isA ;
}

PBWT *pbwtRead (FILE *fp) 
{
  int m, n ;
  PBWT *p ;
  static char tag[5] = "test" ;

  if (fread (tag, 1, 4, fp) != 4 || (strcmp (tag, "PBWT") && strcmp (tag, "GBWT"))) /* early versions wrote GBWT */
    die ("failed to recognise file type in pbwtRead - was it written by pbwt?") ;
  if (fread (&m, sizeof(int), 1, fp) != 1)
    die ("error reading n in pbwtRead") ;
  if (fread (&n, sizeof(int), 1, fp) != 1)
    die ("error reading n in pbwtRead") ;
  p = pbwtCreate (m) ;
  p->N = n ;
  if (fread (&n, sizeof(int), 1, fp) != 1)
    die ("error reading pbwt file") ;
  p->yz = arrayCreate (n, uchar) ;
  array(p->yz, n-1, uchar) = 0 ; /* sets arrayMax */
  if (fread (arrp(p->yz, 0, uchar), sizeof(uchar), n, fp) != n)
    die ("error reading data in pbwt file") ;

  fprintf (stderr, "read pbwt file %d bytes: M, N are %d, %d\n", n, p->M, p->N) ;
  return p ;
}

static BOOL readMatchChrom (PBWT *p, FILE *fp)
{
  char *chrom = fgetword (fp) ;

  if (strcmp (chrom, "."))
    if (!p->chrom) 
      p->chrom = strdup (chrom) ;
    else if (strcmp (chrom, p->chrom))
      return FALSE ;
  return TRUE ;
}

void pbwtReadSites (PBWT *p, FILE *fp)
{
  int n = 0, i ;
  char *chrom, c ;
  Site *s ;
  int line = 1 ;
  Array varTextArray = arrayCreate (256, char) ;

  if (!p) die ("pbwtReadSites called without a valid pbwt") ;

  p->sites = arrayReCreate (p->sites, p->N, Site) ;
  while (!feof(fp))
    if (readMatchChrom (p, fp))	/* if p->chrom then match, else set if not "." */
      { if (feof(fp)) break ;
	s = arrayp(p->sites, n++, Site) ;
	s->x = 0 ; while (isdigit(c = fgetc(fp))) s->x = s->x*10 + (c-'0') ;
	if (!feof(fp) && c != '\n')
	  { if (!isspace (c)) die ("bad position line %d in sites file", line) ;
	    while (isspace(c = fgetc(fp)) && c != '\n') ;
	    if (c == '\n') die ("bad end of line at line %d in sites file", line) ;
	    i = 0 ; array(varTextArray, i++, char) = c ;
	    while ((c = fgetc(fp)) && c != '\n') array(varTextArray, i++, char) = c ;
	    array(varTextArray, i, char) = 0 ;
	    if (!p->variationDict) p->variationDict = dictCreate(32) ;
	    dictAdd (p->variationDict, arrayp(varTextArray,0,char), &s->varD) ;
	  }
	++line ;
      }
    else if (!feof(fp))
      die ("failed to match chromosome in sites file: line %d", line) ;

  if (ferror (fp)) die ("error reading sites file") ;
  
  if (n != p->N)
    { fprintf (stderr, "sites file contains %d sites not %d as in pbwt - reject these sites\n",
	       n, p->N) ;
      arrayDestroy (p->sites) ; p->sites = 0 ;
    }
  else
    fprintf (stderr, "read %d sites from file\n", n) ;

  arrayDestroy (varTextArray) ;
}

void pbwtReadSamples (PBWT *p, FILE *fp) /* for now assume all samples diploid */
{
  char *name, c ;
  int line = 0 ;
  int i = 0, n ;
  Array nameArray = arrayCreate (1024, char) ;

  if (!p) die ("pbwtReadSamples called without a valid pbwt") ;

  p->samples = arrayReCreate (p->samples, p->M/2, Site) ;
  if (!p->sampleNameDict) p->sampleNameDict = dictCreate (p->M/2) ;

  while (!feof(fp))
    { ++line ;
      n = 0 ;
      while ((c = getc(fp)) && !isspace(c) && !feof (fp)) array(nameArray, n++, char) = c ;
      if (feof(fp)) break ;
      if (!n) die ("no name line %d in samples file", line) ;
      array(nameArray, n++, char) = 0 ;
      dictAdd (p->sampleNameDict, arrp(nameArray,0,char), 
	       &(arrayp(p->samples, i, Sample)->nameD)) ;
      arrayp(p->samples, i+1, Sample)->nameD = arrp(p->samples, i, Sample)->nameD ;
      i += 2 ;
      /* now remove the rest of the line, for now */
      while (c != '\n' && !feof(fp)) c = getc(fp) ;
    }

  fprintf (stderr, "read %d sample names\n", dictMax(p->sampleNameDict)) ;

  if (i > p->M) die ("too many samples %d > %d in readSamples - for now assume diploid", i/2, p->M/2) ;

  arrayDestroy (nameArray) ;
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
  uchar *x, *yz, *y2 ;		/* original, sorted, compressed, uncompressed */
  int *a ;
  Update *u ;

  parseMacsHeader (fp, &M, &L) ;
  p = pbwtCreate (M) ;
  u = updateCreate (M, 0) ;
  x = myalloc (M+1, uchar) ; x[M] = Y_SENTINEL ;	/* sentinel required for packing */
  yz = myalloc (M, uchar) ;
  y2 = myalloc (M, uchar) ;

  p->sites = arrayCreate(4096, Site) ;
  p->yz = arrayCreate(4096*32, uchar) ;
  p->N = 0 ;
  while (parseMacsSite (fp, arrayp(p->sites,p->N,Site), M, L, x))
    { for (j = 0 ; j < M ; ++j) u->y[j] = x[u->a[j]] ; /* next character in sort order: BWT */
      nyPack = pack3 (u->y, M, yz) ;
      for (j = 0 ; j < nyPack ; ++j) array(p->yz,arrayMax(p->yz),uchar) = yz[j] ;
      updateForwardsA (u) ;
      if (isCheck)
	{ int nx0 = 0 ;
	  nUnpack = unpack3 (yz, M, y2, &n0) ;
	  for (j = 0 ; j < M ; j++)
	    if (y2[j] != u->y[j])
	      die ("Mismatch y2:%d != y:%d at %d line %d", y2[j], u->y[j], j, i) ;
	  if (nUnpack != nyPack)
	    die ("Unpack %d is not equal to pack %d line %d", nUnpack, nyPack, i) ;
	  for (j = 0 ; j < M ; ++j) if (!x[j]) ++nx0 ;
	  if (nx0 != n0)
	    die ("Mismatch nx0:%d != n0:%d line %d", nx0, n0) ;
      	}
      if (isStats)
	{ if (!isCheck) nUnpack = unpack3 (yz, M, y2, &n0) ;
	  nxPack = pack3 (x, M, yz) ;
	  printf ("Site\t%d\t%d\t%d\t%d\t%d\n", 
		  p->N, arrp(p->sites,p->N,Site)->x, M-n0, nxPack, nyPack) ;
	  nxTot += nxPack ;
	  nyTot += nyPack ;
      	}
      p->N++ ;
      if (nCheckPoint && !(p->N % nCheckPoint))
	pbwtCheckPoint (p) ;
    }

  fprintf (stderr, "read MaCS file: M, N are\t%d\t%d\n", M, p->N) ;
  if (isStats)
    fprintf (stderr, "                xtot, ytot are\t%d\t%d\n", nxTot, nyTot) ;

  free(x) ; free (yz) ; free (y2) ; updateDestroy (u) ;

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

static BOOL parseVcfqSite (PBWT **pp, FILE *fp, Array x) /* parse VCFQ line */
{
  char c, *chrom, *var ;
  int pos, m, len ;
  Site *s ;
  PBWT *p = *pp ;

  if (!p)			/* first call on a file - don't know M yet */
    chrom = strdup (fgetword (fp)) ;
  else
    while (!feof (fp) && !readMatchChrom (p, fp))
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
    { p = *pp = pbwtCreate (m) ;
      if (strcmp (chrom, ".")) p->chrom = chrom ;
      p->variationDict = dictCreate (32) ;
      p->sites = arrayCreate(4096, Site) ;
      array(x,p->M,uchar) = Y_SENTINEL ; /* sentinel required for packing */
    }
  s = arrayp(p->sites, arrayMax(p->sites), Site) ;
  s->x = pos ;
  dictAdd (p->variationDict, var, &s->varD) ;

  ++p->N ;
  return TRUE ;
}

PBWT *pbwtReadVcfq (FILE *fp)
{
  PBWT *p = 0 ;
  int j, nyPack ;
  uchar *x, *yz ;		/* original, sorted, compressed */
  int *a ;
  Array xArray = arrayCreate (10000, uchar) ;
  Update *u ;

  while (parseVcfqSite (&p, fp, xArray)) /* create p first time round */
    { if (!p->yz)		/* first line; p was just made! */
	{ p->yz = arrayCreate(4096*32, uchar) ;
	  yz = myalloc (p->M, uchar) ;
	  u = updateCreate (p->M, 0) ;
	}
      x = arrp(xArray,0,uchar) ;
      for (j = 0 ; j < p->M ; ++j) u->y[j] = x[u->a[j]] ;
      nyPack = pack3 (u->y, p->M, yz) ;
      for (j = 0 ; j < nyPack ; ++j) array(p->yz,arrayMax(p->yz),uchar) = yz[j] ;
      updateForwardsA (u) ;
      if (nCheckPoint && !(p->N % nCheckPoint))
	pbwtCheckPoint (p) ;
    }

  fprintf (stderr, "read vcfq file") ;
  if (p->chrom) fprintf (stderr, " for chromosome %s", p->chrom) ;
  fprintf (stderr, ": M, N are\t%d\t%d; yz length is %d\n", p->M, p->N, arrayMax(p->yz)) ;

  arrayDestroy(xArray) ; free (yz) ; updateDestroy (u) ;

  return p ;
}

/*************** write haplotypes ******************/

void pbwtWriteHaplotypes (FILE *fp, PBWT *p)
{
  int i, j, n = 0, M = p->M ;
  uchar *hap = myalloc (M, uchar) ;
  Update *u = updateCreate (M, 0) ;

  for (i = 0 ; i < p->N ; ++i)
    { n += unpack3 (arrp(p->yz,n,uchar), M, u->y, 0) ;
      for (j = 0 ; j < M ; ++j) hap[u->a[j]] = u->y[j] ;
      for (j = 0 ; j < M ; ++j) putc (hap[j]?'1':'0', fp) ;
      putc ('\n', fp) ; fflush (fp) ;
      updateForwardsA (u) ;
    }
  free (hap) ; updateDestroy (u) ;

  fprintf (stderr, "written haplotype file: %d rows of %d\n", p->N, M) ;
}

/******************* end of file *******************/

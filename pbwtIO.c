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
 * Last edited: Mar 11 16:38 2016 (rd)
 * readPhase updated for chromopainter and chromopainter v2 formats
 * Created: Thu Apr  4 11:42:08 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"
#include <ctype.h>

int nCheckPoint = 0 ;	/* if set non-zero write pbwt and sites files every n sites when parsing external files */

static BOOL isWriteImputeRef = FALSE ;	/* modifies WriteSites() and WriteHaplotypes() for pbwtWriteImputeRef */

/* basic function to store packed PBWT */

void pbwtWrite (PBWT *p, FILE *fp) /* just writes compressed pbwt in yz */
{
  if (!p || !p->yz) die ("pbwtWrite called without a valid pbwt") ;
  if (!p->aFstart || !p->aFend) die ("pbwtWrite called without start and end indexes") ;
  /* version 2 added start and end indexes */
  if (fwrite ("PBW3", 1, 4, fp) != 4) /* version 3 with 8 byte pbwt size */
    die ("error writing PBWT in pbwtWrite") ;
  if (fwrite (&p->M, sizeof(int), 1, fp) != 1)
    die ("error writing M in pbwtWrite") ;
  if (fwrite (&p->N, sizeof(int), 1, fp) != 1)
    die ("error writing N in pbwtWrite") ;
  if (fwrite (p->aFstart, sizeof(int), p->M, fp) != p->M)
    die ("error writing aFstart in pbwtWrite") ;
  if (fwrite (p->aFend, sizeof(int), p->M, fp) != p->M)
    die ("error writing aFend in pbwtWrite") ;
  long n = arrayMax(p->yz) ;
  if (fwrite (&n, sizeof(long), 1, fp) != 1)
    die ("error writing n in pbwtWrite") ;
  if (fwrite ("    ", 1, 4, fp) != 4)
    die ("error writing padding space in pbwtWrite") ;
  if (fwrite (arrp(p->yz, 0, uchar), sizeof(uchar), arrayMax(p->yz), fp) != arrayMax(p->yz))
    die ("error writing data in pbwtWrite") ;

  fprintf (logFile, "written %ld chars pbwt: M, N are %d, %d\n", arrayMax(p->yz), p->M, p->N) ;
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

  fprintf (logFile, "written %d sites from %d to %d\n", p->N, 
	   arrp(p->sites, 0, Site)->x, arrp(p->sites, p->N-1, Site)->x) ;
}

void pbwtWriteSamples (PBWT *p, FILE *fp)
{
  if (!p || !p->samples) die ("pbwtWriteSamples called without samples") ;

  static const char *typeStr = "\0ifZ" ;
  int i, j, count = 0 ;
  for (i = 0 ; i < p->M ; i += pbwtSamplePloidy (p, i), count++)
    { Sample *s = sample (arr(p->samples, i, int)) ;
      writeMetaData(sampleName(s), s->metaData, fp) ;
    }     
  if (ferror (fp)) die ("error writing samples file") ;

  fprintf (logFile, "written %d samples\n", count) ;
}

void writeDataOffset (FILE *fp, char *name, Array data, Array offset, int N)
{
  if (!offset || !data) die ("write %s called without data", name) ;
  int dummy = -1 ;	/* ugly hack to mark that we now use longs not ints */
  if (fwrite (&dummy, sizeof(int), 1, fp) != 1)
    die ("error writing marker in write %s", name) ;
  long n = arrayMax(data) ;
  if (fwrite (&n, sizeof(long), 1, fp) != 1)
    die ("error writing n in write %s", name) ;
  if (fwrite (arrp(data, 0, uchar), sizeof(uchar), n, fp) != n)
    die ("error writing data in write %s", name) ;
  if (fwrite (arrp(offset, 0, long), sizeof(long), N, fp) != N)
    die ("error writing offsets in write %s", name) ;

  fprintf (logFile, "written %ld chars compressed %s data\n", n, name) ;
}

void pbwtWriteMissing (PBWT *p, FILE *fp)
{ writeDataOffset (fp, "missing", p->zMissing, p->missingOffset, p->N) ; }

void pbwtWriteDosage (PBWT *p, FILE *fp)
{ writeDataOffset (fp, "dosage", p->zDosage, p->dosageOffset, p->N) ; }

void pbwtWriteReverse (PBWT *p, FILE *fp)
{
  if (!p || !p->zz) die ("pbwtWriteReverse called without reverse pbwt") ;

  Array tz = p->yz ; p->yz = p->zz ;
  int* tstart = p->aFstart ; p->aFstart = p->aRstart ;
  int* tend = p->aFend ; p->aFend = p->aRend ;

  fprintf (logFile, "reverse: ") ; pbwtWrite (p, fp) ;
  
  p->yz = tz ; p->aFstart = tstart ; p->aFend = tend ;
}


void pbwtWriteAll (PBWT *p, char *root)
{
  FILE *fp ;
#define FOPEN_W(tag)  if (!(fp = fopenTag (root, tag, "w"))) die ("failed to open root.%s", tag)
  FOPEN_W("pbwt") ; pbwtWrite (p, fp) ; fclose (fp) ;
  if (p->sites) { FOPEN_W("sites") ; pbwtWriteSites (p, fp) ; fclose (fp) ; }
  if (p->samples) { FOPEN_W("samples") ; pbwtWriteSamples (p, fp) ; fclose (fp) ; }
  if (p->missingOffset) { FOPEN_W("missing") ; pbwtWriteMissing (p, fp) ; fclose (fp) ; }
  if (p->dosageOffset) { FOPEN_W("dosage") ; pbwtWriteDosage (p, fp) ; fclose (fp) ; }
  if (p->zz) { FOPEN_W("reverse") ; pbwtWriteReverse (p, fp) ; fclose (fp) ; }
}

void pbwtWritePhase (PBWT *p, char *filename)
{
  FILE *fp ;
  if (!(fp = fopen (filename, "w"))) die ("failed to open %s",filename);
  fprintf(fp,"%i\n",p->M);
  fprintf(fp,"%i\nP",p->N);
  int i ; for (i = 0 ; i < p->N ; i++) fprintf(fp," %i",arrayp(p->sites,i,Site)->x);
  fprintf(fp,"\n");

  pbwtWriteTransposedHaplotypes(p,fp);
}

void pbwtCheckPoint (PbwtCursor *u, PBWT *p)
{
  static BOOL isA = TRUE ;
  char fileNameRoot[20] ;

  pbwtCursorToAFend (u, p) ;
  sprintf (fileNameRoot, "check_%c", isA ? 'A' : 'B') ;
  pbwtWriteAll (p, fileNameRoot) ;

  isA = !isA ;
}

/*******************************/

PBWT *pbwtRead (FILE *fp) 
{
  int m, n ;
  long nz ;
  PBWT *p ;
  static char tag[5] = "test" ;
  char pad[4] ;
  int version ;

  if (fread (tag, 1, 4, fp) != 4) die ("failed to read 4 char tag - is file readable?") ;
  if (!strcmp (tag, "PBW3")) version = 3 ; /* current version */
  else if (!strcmp (tag, "PBW2")) version = 2 ; /* with 4 byte count */
  else if (!strcmp (tag, "PBWT")) version = 1 ; /* without start, end indexes */
  else if (!strcmp (tag, "GBWT")) version = 0 ; /* earliest version */
  else die ("failed to recognise file type %s in pbwtRead - was it written by pbwt?", tag) ;

  if (fread (&m, sizeof(int), 1, fp) != 1) die ("error reading m in pbwtRead") ;
  if (fread (&n, sizeof(int), 1, fp) != 1) die ("error reading n in pbwtRead") ;
  p = pbwtCreate (m, n) ;
  free(p->aFstart);
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

  if (version <= 2)
    { if (fread (&n, sizeof(int), 1, fp) != 1) die ("error reading pbwt file") ;
      nz = n ;
    }
  else
    if (fread (&nz, sizeof(long), 1, fp) != 1 ||
	fread (pad, 1, 4, fp) != 4) die ("error reading pbwt file") ;

  p->yz = arrayCreate (nz, uchar) ;
  array(p->yz, nz-1, uchar) = 0 ; /* sets arrayMax */
  if (fread (arrp(p->yz, 0, uchar), sizeof(uchar), nz, fp) != nz)
    die ("error reading data in pbwt file") ;

  fprintf (logFile, "read pbwt %s file with %ld bytes: M, N are %d, %d\n", tag, nz, p->M, p->N) ;
  return p ;
}

static BOOL readMatchChrom (char **pChrom, FILE *fp)
{
  char *newChrom = fgetword (fp) ;

  if (strcmp (newChrom, "."))
    { if (!*pChrom) 
	*pChrom = strdup (newChrom) ;
      else if (strcmp (newChrom, *pChrom)) 
	return FALSE ;
    }
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
  
  fprintf (logFile, "read %ld sites on chromosome %s from file\n", arrayMax(sites), *chrom) ;

  arrayDestroy (varTextArray) ;
  return sites ;
}

void pbwtReadSites (PBWT *p, FILE *fp)
{
  if (!p) die ("pbwtReadSites called without a valid pbwt") ;

  p->sites = pbwtReadSitesFile (fp, &p->chrom) ;
  if (arrayMax(p->sites) != p->N)
    die ("sites file contains %ld sites not %d as in pbwt", arrayMax(p->sites), p->N) ;
}

Array readRefFreq (PBWT *p, FILE *fp)
{
  Array a = arrayCreate (p->N, Site) ;
  char chrom[256], var[256] ;
  int pos ;
  double freq ;
  while (!feof(fp))
    { if (fscanf (fp, "%s\t%d\t%lf\t", chrom, &pos, &freq) != 3 && !feof (fp))
	die ("can't read ref freq file properly") ;
      if (strcmp (chrom, p->chrom))
	die ("chromosome mismatch in readRefFreq '%s' is not '%s'", chrom, p->chrom) ;
      char *cp = var ; while ((*cp = getc(fp)) != '\n' && !feof(fp)) ++cp ; *cp = 0 ;
      Site *s = arrayp(a, arrayMax(a), Site) ;
      s->x = pos ; 
      s->refFreq = freq ;
      dictAdd (variationDict, var, &s->varD) ;
    }
  return a ;
}

void pbwtReadRefFreq (PBWT *p, FILE *fp)
{ if (!p || !p->sites) die ("pbwtReadRefFreq called without current site information") ;
  Array a = readRefFreq (p, fp) ;
  
  int i = 0, j = 0 ;
  Site *ps = arrp(p->sites,i,Site), *as = arrp(a,j,Site) ;
  while (i < p->N)
    { while (j < arrayMax(a) && 
	     (as->x < ps->x || (as->x == ps->x && as->varD < ps->varD))) { ++j ; ++as ; }
      if (ps->x == as->x && ps->varD == as->varD) ps->refFreq = as->refFreq ;
      ++i ; ++ps ;
    }
}

// Array pbwtReadSamplesFile (FILE *fp) /* for now assume all samples diploid */
// /* should add code to this to read father and mother and population
//    propose to use IMPUTE2 format for this */
// {
//   char *name, c ;
//   int line = 0 ;
//   Array nameArray = arrayCreate (1024, char) ;
//   Array samples = arrayCreate (1024, int) ;

//   while (!feof(fp))
//     { int n = 0 ;
//       while ((c = getc(fp)) && !isspace(c) && !feof (fp)) array(nameArray, n++, char) = c ;
//       if (feof(fp)) break ;
//       if (!n) die ("no name line %ld in samples file", arrayMax(samples)+1) ;
//       array(nameArray, n++, char) = 0 ;
//      /* next bit of code deals with header lines for IMPUTE2 samples file, so can read that */
//       if (!strcmp (arrp(nameArray,0,char), "ID_1") && !arrayMax(samples))
// 	{ while ((c = getc(fp)) && c != '\n' && !feof(fp)) ; /* remove header line */
// 	  while ((c = getc(fp)) && c != '\n' && !feof(fp)) ; /* and next line of zeroes?? */
// 	  continue ;
// 	}
//       array(samples,arrayMax(samples),int) = sampleAdd(arrp(nameArray,0,char), 0, 0, 0, 0, 0)->nameD ;
//       /* now remove the rest of the line, for now */
//       while (c != '\n' && !feof(fp)) c = getc(fp) ;
//     }
//   arrayDestroy (nameArray) ;

//   fprintf (logFile, "read %ld sample names\n", arrayMax(samples)) ;

//   return samples ;
// }

  /*
     FMF format - tab delimited flat metadata format where first column is the 
     sample name followed by metadata in the format key:type:value with type 
     being Z for a string, f for a real number and i for an integer:
      sample1   sex:Z:M    height:f:1.73     region:Z:WestEurasia     foo:i:10
      sample2   sex:Z:F    height:f:1.64     region:Z:WestEurasia     bar:i:20
  */

  /*
     FAM format - first 6 columns of the PED format (used for binary PED format). space delimited with no header:
     1. Family ID
     2. Individual ID
     3. Paternal ID
     4. Maternal ID
     5. Sex (1=male; 2=female; other=unknown)
     6. Phenotype (-1=missing; 0=missing; 1=unaffected; 2=affected)
  */

  /*
     HAPS/SAMPLE and GEN/SAMPLE format - space delimited with at least 3 columns with 2 header lines:
     1. individual ID 1 (ID_1)
     2. individual ID 2 (ID_2)
     3. missing data proportion (missing)
     4-. ignored, covariates named in the first header line and type defined in the second header line
         (0=for first 3 columns; D=discrete covariate (unsigned int); C=continuous covariate (float); P=continuous phenotype (float); B=binary phenotype (0 or 1))
  */

  /*
     HAPS/LEGEND/SAMPLE format - space delimited with at least 4 columns with 1 header line :
     1. individual ID (sample)
     2. group ID 1 (population)
     3. group ID 2 (group)
     4. sex (1=male; 2=female; other=unknown)
     5-. ignored
  */
     
Array pbwtReadSamplesFile (FILE *fp) {
  // populate sample fields from FMF sample file
  char *line = calloc(1024, sizeof(char)); line[0] = '\0';
  char *key = calloc(128, sizeof(char)); key[0] = '\0';
  char *value = calloc(128, sizeof(char)); value[0] = '\0';
  char *name = calloc(128, sizeof(char)); name[0] = '\0';
  char *str;
  char type;

  Array samples = arrayCreate (1024, int) ;

  while (1) {
    char *rmme = line;
    line = fgets(line, 1024, fp);
    if ( !line ) free(rmme);
    if(feof(fp)) break;
    else if(sscanf(line, "ID_1%s", name)) { // IMPUTE2 file
      die("IMPUTE2 style sample files not yet supported - coming soon");
    }
    else { // FMF file (single column valid)
      str = strtok(line, "\t\n"); // test it works with single column
      strcpy(name, str);
      int v, k = array(samples,arrayMax(samples),int) = sampleAdd(name) ;
      Sample *s = sample (k) ;
      while(str = strtok(0, "\t\n"))
      {
        if(sscanf(str, "%127[^:]:%1[^:]:%127[^:]", key, &type, value) == 3)
        {
            addMetaData(s->metaData, key, value, type) ;
        }
        if(!strcasecmp(key, "mother")) { mother (s, value) ; }
        if(!strcasecmp(key, "father")) { father (s, value) ; }
        if(!strcasecmp(key, "family")||!strcasecmp(key, "familyID")) { familyName (s, value) ; }
        if(!strcasecmp(key, "pop")||!strcasecmp(key, "population")) { popName (s, value) ; }
        if(!strcasecmp(key, "gender")||!strcasecmp(key, "sex"))
        {
          if(!strcasecmp(value, "M")||!strcasecmp(value, "male")) s->isMale = TRUE ;
          if(!strcasecmp(value, "F")||!strcasecmp(value, "female")) s->isFemale = TRUE ;
        }
      }
    }    
  }
  free(line); free(key); free(value); free(name);
  fprintf (logFile, "read %ld sample names\n", arrayMax(samples)) ;
  return samples ;
}

void pbwtReadSamples (PBWT *p, FILE *fp)
{
  if (!p) die ("pbwtReadSamples called without a valid pbwt") ;

  // NB: pbwtReadSamplesFile() calls sampleAdd() which uses static sample
  // dictionary and therefore the sample indexes will be offset when pbwtReadSamples()
  // is called multiple times.
  int ioff = sampleCount();

  Array samples = pbwtReadSamplesFile (fp) ;
  int i,j, count_x2 = 0, count_x1 = 0, count_y1 = 0, count_y0 = 0, count2 = 0; 
  p->samples = arrayReCreate(p->samples, p->M, int) ;

  for (i=0; i<arrayMax(samples); i++)
  {
      Sample *s = sample (i+ioff) ;
      if ( s->isFemale ) { count_x2++; count_y0++; count2++; }
      else { count_x1++; count_y1++; count2++; }
  }

  if ( 2*count2 == p->M )
      fprintf (logFile, "Note: number of haplotypes (%ld) is consistent with %d diploid samples\n", arrayMax(samples), p->M) ;
  else if ( 2*count_x2 + count_x1 == p->M )
  {
      p->isX = TRUE;
      fprintf (logFile, "Note: number of haplotypes (%d) is consistent with %d haploid and %d diploid samples\n", p->M, count_x1,count_x2) ;
  }
  else if ( count_y1 == p->M )
  {
      p->isY = TRUE;
      fprintf (logFile, "Note: number of haplotypes (%d) consistent with %d haploid samples and %d samples with zero ploidy\n", p->M, count_y1,count_y0) ;
  }
  else
      die ("number of haplotypes (%d) is not consistent with the number of samples (%ld): expected %d, %d or %d haplotypes", p->M, arrayMax(samples),count2,2*count_x2 + count_x1,count_y1) ;

  for (i=0,j=0; i<arrayMax(samples); i++)
  {
      Sample *s = sample (i+ioff) ;
      if (p->isY && s->isFemale) continue ;
      array(p->samples, j++, int) = arr(samples, i, int) ;
      if (p->isX && s->isMale) continue ;
      array(p->samples, j++, int) = arr(samples, i, int) ;
  }
  arrayDestroy (samples) ;
}

static void readDataOffset (FILE *fp, char *name, Array *data, Array *offset, int N)
{
  long n ;			/* size of data file */
  int dummy ; 
  if (fread (&dummy, sizeof(int), 1, fp) != 1) 
    die ("read error in read %s", name) ;
  if (dummy != -1) n = dummy ;	/* old version with ints not longs */
  else if (fread (&n, sizeof(long), 1, fp) != 1) 
    die ("read error in read %s", name) ;

  *data = arrayReCreate (*data, n, uchar) ;
  if (fread (arrp(*data, 0, uchar), sizeof(uchar), n, fp) != n)
    die ("error reading z%s in pbwtRead%s", name, name) ;
  arrayMax(*data) = n ;
  fprintf (logFile, "read %ld chars compressed %s data\n", n, name) ;

  *offset = arrayReCreate (*offset, N, long) ;
  if (dummy != -1)		/* old version with ints not longs */
    { /* abuse *offset to hold ints, then update in place */
      if (fread (arrp(*offset, 0, int), sizeof(int), N, fp) != N) 
	die ("error reading %s in pbwtRead%s", name, name) ;
      for (n = N ; n-- ;) 
	arr(*offset, n, long) = arr(*offset, n, int) ; /* !! */
    }
  else
    if (fread (arrp(*offset, 0, long), sizeof(long), N, fp) != N)
      die ("error reading %s in pbwtRead%s", name, name) ;
  arrayMax(*offset) = N ;
}

void pbwtReadMissing (PBWT *p, FILE *fp)
{ readDataOffset (fp, "missing", &p->zMissing, &p->missingOffset, p->N) ; }

void pbwtReadDosage (PBWT *p, FILE *fp)
{ readDataOffset (fp, "dosage", &p->zDosage, &p->dosageOffset, p->N) ; }

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

PBWT *pbwtReadAll (char *root)
{
  PBWT *p ;

  FILE *fp ;
  if ((fp = fopenTag (root, "pbwt", "r"))) { p = pbwtRead (fp) ; fclose (fp) ; } 
  else die ("failed to open %s.pbwt", root) ;
  if ((fp = fopenTag (root, "sites","r")))  { pbwtReadSites (p, fp) ; fclose (fp) ; }
  if ((fp = fopenTag (root, "samples","r"))) { pbwtReadSamples (p, fp) ; fclose (fp) ; }
  if ((fp = fopenTag (root, "missing","r"))) { pbwtReadMissing (p, fp) ; fclose (fp) ; }
  if ((fp = fopenTag (root, "dosage","r"))) { pbwtReadDosage (p, fp) ; fclose (fp) ; }
  if ((fp = fopenTag (root, "reverse","r"))) { pbwtReadReverse (p, fp) ; fclose (fp) ; }

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
  double dummy = atof(fgetword(fp)) ;		/* ignore the time - dummy stops compile warning */
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
      if (nCheckPoint && !(p->N % nCheckPoint)) pbwtCheckPoint (u, p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (logFile, "read MaCS file: M, N are\t%d\t%d\n", M, p->N) ;
  if (isStats)
    fprintf (logFile, "                xtot, ytot are\t%d\t%d\n", nxTot, nyTot) ;

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
      if (nCheckPoint && !(p->N % nCheckPoint))	pbwtCheckPoint (u, p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (logFile, "read %s file", type) ;
  if (p->chrom) fprintf (logFile, " for chromosome %s", p->chrom) ;
  fprintf (logFile, ": M, N are\t%d\t%d; yz length is %ld\n", p->M, p->N, arrayMax(p->yz)) ;

  arrayDestroy(xArray) ; pbwtCursorDestroy (u) ;

  return p ;
}

typedef BOOL (*ParseTwoLineFunc)(PBWT** pp, FILE *f, FILE *l, Array a) ;

static PBWT *pbwtReadLineTwoFile (FILE *fp, FILE *lp, char* type, ParseTwoLineFunc parseLine)
{
  PBWT *p = 0 ;
  int j ;
  uchar *x ;		/* original, sorted, compressed */
  int *a ;
  Array xArray = arrayCreate (10000, uchar) ;
  PbwtCursor *u ;

  /* skip header of legend file */
  while (!feof(lp))
  {  char c = getc(lp) ; if (c == '\n') break ;
     }

  while ((*parseLine) (&p, fp, lp, xArray)) /* create p first time round */
    { if (!p->yz)		/* first line; p was just made! */
	{ p->yz = arrayCreate(4096*32, uchar) ;
	  u = pbwtCursorCreate (p, TRUE, TRUE) ;
	}
      x = arrp(xArray,0,uchar) ;
      for (j = 0 ; j < p->M ; ++j) u->y[j] = x[u->a[j]] ;
      pbwtCursorWriteForwards (u) ;
      if (nCheckPoint && !(p->N % nCheckPoint))	pbwtCheckPoint (u, p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (stderr, "read %s file", type) ;
  if (p->chrom) fprintf (stderr, " for chromosome %s", p->chrom) ;
  fprintf (stderr, ": M, N are\t%d\t%d; yz length is %ld\n", p->M, p->N, arrayMax(p->yz)) ;

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
  
  fgetword(fp) ;	fgetword(fp) ;	/* ignore first two name fields */
  
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


static BOOL parseHapLegendLine (PBWT **pp, FILE *fp, FILE *lp, Array x) /* same as parseGenLine - slightly simpler */
{
  PBWT *p = *pp ;
  
  fgetword(lp) ;	/* ignore first name field */
  
  int pos = atoi (fgetword (lp)) ;
  char *var = getVariation (lp) ; /* but need to change ' ' separator to '\t' */
  if (feof (lp)) return FALSE ;
  char *cp = var ; while (*cp && *cp != ' ') ++cp ; if (*cp == ' ') *cp = '\t' ; else die ("missing separator in line %d, var is %d", p?p->N:0, var) ;
  while (!feof(lp))
  {  char c = getc(lp) ; if (c == '\n') break ;
     }

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
  if (nGenMissing) fprintf (logFile, "%ld missing genotypes set to 00\n", nGenMissing) ;
  return p ;
}

PBWT *pbwtReadHap (FILE *fp, char *chrom)
{
  PBWT *p = pbwtReadLineFile (fp, "hap", parseHapLine) ;
  p->chrom = strdup (chrom) ;
  return p ;
}

PBWT *pbwtReadHapLegend (FILE *fp, FILE *lp, char *chrom)
{
  PBWT *p = pbwtReadLineTwoFile (fp, lp, "hap-legend", parseHapLegendLine) ;
  p->chrom = strdup (chrom) ;
  return p ;
}

PBWT *pbwtReadPhase (FILE *fp) /* Li and Stephens PHASE format */
{
  int nhaps=0,nsnps=0,ninds=0;
  int version=2;
  int l1=atoi(fgetword (fp)) ; if (getc(fp) != '\n') die ("bad first line in phase file") ;
  int l2 = atoi (fgetword (fp)) ;  if (getc(fp) != '\n') die ("bad second line in phase file") ;
  char * tc=fgetword (fp);
  int l3 = atoi (tc) ;  
  if (tc[0] == 'P') { // version 2
    nhaps=l1;
    nsnps=l2;
    ninds=nhaps/2;
  }else{ // version 1
    if (getc(fp) != '\n') die ("bad third line in phase file") ;
    ninds=l2;
    nhaps=l2*2;
    nsnps=l3;
    fgetword (fp); /* Remove the initial P */
    version=1;
  }
  fprintf(logFile,"Reading %i SNPs %i haplotypes and %i individuals from PHASE format version %i\n",nsnps,nhaps,ninds,version);
  PBWT *p = pbwtCreate (nhaps, nsnps) ;
  p->chrom = strdup ("0") ; /* chromosome number not available from phase files */
  p->sites = arrayCreate (4096, Site) ;
  int i ; for (i = 0 ; i < p->N ; ++i) {
    tc=fgetword(fp);
    arrayp(p->sites,i,Site)->x = atoi (tc) ;
  }
  if (getc(fp) != '\n') die ("bad location line in phase file") ;
  char var[2] ="S\0"; var[1] = 0 ; /*dummy variation line - used to create variationDict */
  for (i = 0 ; i < p->N ; ++i) 
    {
      if(version==1) *var = getc(fp) ;
      dictAdd (variationDict, var, &(arrayp(p->sites,i,Site)->varD)) ; }
  if(version==1) if (getc(fp) != '\n') die ("bad 5th line in phase file") ;

  uchar **data = myalloc (p->N, uchar*) ;
  for (i = 0 ; i < p->N ; ++i) data[i] = myalloc (p->M, uchar) ;
  int j ; for (j = 0 ; j < p->M ; ++j)
    { for (i = 0 ; i < p->N ; ++i) data[i][j] = getc(fp) - '0' ;
      if (getc(fp) != '\n') die ("bad %dth line in phase file", 7+j-version) ;
    }
  p->yz = arrayCreate(4096*32, uchar) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  for (i = 0 ; i < p->N ; ++i)
    { for (j = 0 ; j < p->M ; ++j) u->y[j] = data[i][u->a[j]] ;
      pbwtCursorWriteForwards (u) ;
      if (nCheckPoint && !((i+1) % nCheckPoint)) pbwtCheckPoint (u, p) ;
    }
  pbwtCursorToAFend (u, p) ;

  fprintf (logFile, "read phase file") ;
  if (p->chrom) fprintf (logFile, " for chromosome %s", p->chrom) ;
  fprintf (logFile, ": M, N are\t%d\t%d; yz length is %ld\n", p->M, p->N, arrayMax(p->yz)) ;

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

  fprintf (logFile, "written haplotype file: %d rows of %d\n", p->N, M) ;
}

void pbwtWriteTransposedHaplotypes (PBWT *p, FILE *fp)
{
  int i, j ;
  uchar **hap = pbwtHaplotypes (p) ;
  
  for (j = 0 ; j < p->M ; ++j)
    { for (i = 0 ; i < p->N ; ++i) 
	putc (hap[j][i], fp) ;
      putc ('\n', fp) ; fflush (fp) ;
    }
  
  for (i = 0 ; j < p->M ; ++j) free(hap[j]) ;  free (hap) ;
  
  fprintf (logFile, "written transposed haplotype file: %d rows of %d\n", p->M,p->N) ;
}

/*************** write IMPUTE files ********************/

void pbwtWriteImputeRef (PBWT *p, char *fileNameRoot)
{
  FILE *fp ;

  isWriteImputeRef = TRUE ;

  if (!(fp=fopenTag (fileNameRoot, "imputeHaps","w"))) die ("can't open file") ; 
  pbwtWriteHaplotypes (fp, p) ; fclose (fp) ;

  if (!(fp=fopenTag (fileNameRoot, "imputeLegend","w"))) die ("can't open file") ; 
  fprintf (fp, "rsID\tposition\ta0\ta1\n") ; /* header line */
  pbwtWriteSites (p, fp) ; fclose (fp) ;

  isWriteImputeRef = FALSE ;
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
  BOOL isDosage = p->dosageOffset ? TRUE : FALSE ;
  double *d = 0, *ad = isDosage ? myalloc (p->M, double) : NULL;

  for (i = 0 ; i < p->N ; ++i)
    { Site *s = arrp(p->sites, i, Site) ;
      char *als = strdup( dictName(variationDict, s->varD) ), *ss = als ;
      while ( *ss ) { if ( *ss=='\t' ) *ss = '_' ; ss++ ; }
      fprintf (fp, "%s:%d_%s %s:%d_%s %d", p->chrom, s->x, als, p->chrom, s->x, als, s->x) ;
      ss = als ;
      while ( *ss ) { if ( *ss=='_' ) *ss = ' ' ; ss++ ; }
      fprintf (fp, " %s", als) ;
      free(als) ;
      if (isDosage) d = pbwtDosageRetrieve (p, u, d, i) ;
      
      for (j = 0 ; j < p->M ; ++j)
        { hap[u->a[j]] = u->y[j] ;
          if (isDosage) ad[u->a[j]] = d[j] ;
        }
      if (isDosage)
          for (j = 0 ; j < p->M ; j+=2) 
            fprintf (fp, " %f %f %f", (1-ad[j]) * (1-ad[j+1]), ad[j] + ad[j+1] - 2*ad[j]*ad[j+1], ad[j] * ad[j+1]) ;
      else
          for (j = 0 ; j < p->M ; j+=2) 
            if (hap[j] + hap[j+1] == 0) fprintf (fp, " 1 0 0") ;
            else if (hap[j] + hap[j+1] == 1) fprintf (fp, " 0 1 0") ;
            else fprintf (fp, " 0 0 1") ;
      putc ('\n', fp) ; fflush (fp) ;
      pbwtCursorForwardsRead (u) ;
    }

  if (ad) free(ad) ;
  free (hap) ; pbwtCursorDestroy (u) ;
}

/******************* end of file *******************/

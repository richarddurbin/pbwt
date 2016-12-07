/*  File: utils.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Modified by Daniel Lawson (dan.lawson@bristol.ac.uk) in December 2014, adding gzip output for paintSparse
 *  Copyright (C) Genome Research Limited, 1996-
 *-------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free Software
 * Foundation; either version 2.1 of the License, or (at your option) any later
 * version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * You should have received a copy of the GNU Lesser General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *-------------------------------------------------------------------
 * Description: core utility functions
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 28 14:02 2014 (dl)
 * adding gzip output for paintSparse
 * Created: Thu Aug 15 18:32:26 1996 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include "utils.h"

void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  timeUpdate (stderr) ;

  exit (-1) ;
}

void warn (char *format, ...)
{
  static int count = 0 ; 
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  if (++count > 9) die ("too many errors") ;
}

long int totalAllocated = 0 ;

void *_myalloc (long size)
{
  void *p = (void*) malloc (size) ;
  if (!p) die ("myalloc failure requesting %d bytes", size) ;
  totalAllocated += size ;
  return p ;
}

void *_mycalloc (long number, int size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d of size %d bytes", number, size) ;
  totalAllocated += number*size ;
  return p ;
}

/*************************************************/

FILE *fopenTag (char* root, char* tag, char* mode)
{
  if (strlen (tag) > 30) die ("tag %s in fopenTag too long - should be < 30 chars", tag) ;
  char *fileName = myalloc (strlen (root) + 32, char) ;
  strcpy (fileName, root) ;
  strcat (fileName, ".") ;
  strcat (fileName, tag) ;
  FILE *f = fopen (fileName, mode) ;
  free (fileName) ;
  return f ;
}

gzFile gzopenTag (char* root, char* tag, char* mode)
{
  if (strlen (tag) > 40) die ("tag %s in gzopenTag too long - should be < 30 chars", tag) ;
  char *fileName = myalloc (strlen (root) + 42, char) ;
  strcpy (fileName, root) ;
  strcat (fileName, ".") ;
  strcat (fileName, tag) ;
  gzFile f = gzopen (fileName, mode) ;
  free (fileName) ;
  return f ;
}

/*************************************************/

char *fgetword (FILE *f)	// pass NULL to free alloced memory
{
  int n = 0 ;
  static char *buf = 0 ;
  int bufSize = 64 ;
  char *cp ;

  if (!f) { if (buf) free(buf); buf = NULL; return NULL; }

  if (!buf) buf = myalloc (bufSize, char) ;
  cp = buf ;
  while (!feof (f) && (*cp = getc (f)))
    if (isgraph(*cp) && !isspace(*cp))
      { ++cp ; ++n ;
	if (n >= bufSize)
	  { bufSize *= 2 ;
	    if (!(buf = (char*) realloc (buf, bufSize)))
	      die ("fgetword realloc failure requesting %d bytes", bufSize) ;
	  }
      }
    else
      { while ((isspace(*cp) || !isgraph(*cp)) && *cp != '\n' && !feof(f)) *cp = getc (f) ;
	/* previous line was
	while ((*cp = getc(f)) && (isspace(*cp) || !isgraph(*cp)) && *cp != '\n' && !feof(f)) ;
	 */
	ungetc (*cp, f) ;
	break ;
      }
  *cp = 0 ;
  return buf ;
}

/***************** rusage for timing information ******************/

#include <sys/resource.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

#ifdef RUSAGE_STRUCTURE_DEFINITIONS

struct rusage {
  struct timeval ru_utime; /* user time used */
  struct timeval ru_stime; /* system time used */
  long ru_maxrss;          /* integral max resident set size */
  long ru_ixrss;           /* integral shared text memory size */
  long ru_idrss;           /* integral unshared data size */
  long ru_isrss;           /* integral unshared stack size */
  long ru_minflt;          /* page reclaims */
  long ru_majflt;          /* page faults */
  long ru_nswap;           /* swaps */
  long ru_inblock;         /* block input operations */
  long ru_oublock;         /* block output operations */
  long ru_msgsnd;          /* messages sent */
  long ru_msgrcv;          /* messages received */
  long ru_nsignals;        /* signals received */
  long ru_nvcsw;           /* voluntary context switches */
  long ru_nivcsw;          /* involuntary context switches */
};

struct timeval {
  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
  suseconds_t  tv_usec;  /* and microseconds */
} ;

#endif /* RUSAGE STRUCTURE_DEFINITIONS */

void timeUpdate (FILE *f)
{
  static BOOL isFirst = TRUE ;
  static struct rusage rOld ;
  struct rusage rNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fprintf (f, "\tMemory\t%li", totalAllocated) ;   
      fputc ('\n', f) ;
    }
  else
    isFirst = FALSE ;

  rOld = rNew ;
}

/***************** fmf metadata ******************/

static DICT *keyDict ;
static DICT *valueDict ;

void metaDataInit (void)
{
  keyDict = dictCreate (4096) ;
  valueDict = dictCreate (4096) ;
}

void metaDataDestroy (void)
{
  if (keyDict) dictDestroy(keyDict);
  if (valueDict) dictDestroy(valueDict);
}

int addMetaData (Array metaData, char *key, char *value, char type)
{
  if (!metaData) {
    metaData = arrayCreate(4096, MetaData) ;
  }
  MetaData *m = arrayp(metaData, arrayMax(metaData), MetaData) ;

  int k, v ;
  dictAdd(keyDict, key, &k) ;

  if (type == 'i') { m->type = FMF_INT, m->value.i = strtol(value, NULL, 0); }
  else if (type == 'f') { m->type = FMF_REAL, m->value.r = strtod(value, NULL) ; }
  else if (type == 'Z')
  {
    dictAdd(valueDict, value, &v) ;
    m->type = FMF_STR ; m->value.s = v ;
  }
  else m->type = FMF_FLAG ;

  return k ;
}

void writeMetaData (char *rowname, Array metaData, FILE *fp)
{
  static const char *typeStr = "\0ifZ" ;
  fputs (rowname, fp) ;
  if (metaData) 
    {
      int i ;
      for (i = 0 ; i < arrayMax(metaData) ; i++)
        {
          MetaData *m = arrp(metaData, i, MetaData) ;
          fputc ('\t', fp) ;
          fputs (metaDataKey (m), fp) ;
          if (m->type != FMF_FLAG)
            {
              fputc (':', fp) ;
              fputc (typeStr[m->type], fp) ;
              fputc (':', fp) ;
              if (m->type == FMF_INT) fprintf (fp, "%lld", (long long)m->value.i) ;
              else if (m->type == FMF_REAL) fprintf (fp, "%g", m->value.r) ;
              else fputs (metaDataStrValue (m), fp);
            }
        }
    }
  fputc ('\n', fp) ;
  if (ferror (fp)) die ("error writing metaData file") ;
}

char *metaDataKey (MetaData *m) { return dictName (keyDict, m->key) ; }
char *metaDataStrValue (MetaData *m) { return dictName (valueDict, m->value.s) ; }

/********************* end of file ***********************/

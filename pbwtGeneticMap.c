/*  File: pbwtGeneticMap.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2014
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
 * Description: genetic map code for pbwt package
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  7 16:25 2015 (rd)
 * Created: Thu Oct 16 22:53:32 2014 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

typedef struct {
  char* chrom ;			/* chromosome name */
  Array x ;			/* of int - base pair coordinates */
  Array g ;			/* of double - genetic map position */
  int x0 ;			/* bounds of map in base pair coordinates */
  Array z ;			/* genetic map position every 100bp from x0 */
} GeneticMap ;

static GeneticMap map ;

/****************************/

static void buildMap (void)
{
  map.x0 = arr(map.x, 0, int) ;
  int n = (arr(map.x, arrayMax(map.x)-1, int) - map.x0) / 100 ;
  map.z = arrayReCreate (map.z, n, double) ;
  int i = 0, *mapx = arrp(map.x,0,int) ; 
  double *mapg = arrp(map.g,0,double) ;
  array(map.z,i,double) = 0.0 ;
  while (i++ < n)
    { int xi = map.x0 + 100*i ;
      while (mapx[1] < xi) { ++mapx ; ++mapg ; }
      array(map.z,i,double) = 
	*mapg + (xi - mapx[0]) * (mapg[1] - mapg[0]) / (mapx[1] - mapx[0]) ;
    }
}

/****************************/

void readGeneticMap (FILE *fp)
{
  if (strcmp (fgetword (fp), "Chromosome") || 
      strcmp (fgetword (fp), "Position(bp)") ||
      strcmp (fgetword (fp), "Rate(cM/Mb)") ||
      strcmp (fgetword (fp), "Map(cM)")) die ("bad first line in readGeneticMap") ;

  map.x = arrayReCreate (map.x, 100000, int) ;
  map.g = arrayReCreate (map.g, 100000, double) ;

  char chrom[256] ;
  int n = 0, x ;
  double rate, y, oldRate ;
  while (!feof (fp))
    { if (fscanf (fp, "%s\t%d\t%lf%lf", chrom, &x, &rate, &y) == 4)
	{ if (n) 
	    { array(map.g, n, double) = arr(map.g, n-1, double) +
		(x - arr(map.x, n-1, int)) * oldRate ;
	    }
	  else 
	    { array(map.g, n, double) = 0.0 ;
	      map.chrom = strdup (chrom) ;
	    }
	  array(map.x, n, int) = x ;
	  oldRate = rate * 0.000001 ;
	  ++n ;
	}
      else if (!feof (fp)) die ("bad line %d in genetic map file", n+2) ;
    }
  if (!n) die ("no data lines in genetic map file") ;
  if (n==1) die ("only one data line in genetic map file") ;

  buildMap () ;

  fprintf (logFile, "read %d genetic map entries from %d, %f to %d, %f\n",
	   n, arr(map.x, 0, int), arr(map.g, 0, double), 
	   arr(map.x, n-1, int), arr(map.g, n-1, double)) ;
}

/*****************************************/

double geneticMap (int x)
{
  x -= map.x0 ;
  if (x <= 0) return 0.0 ;
  int xi = x / 100 ;
  if (xi >= arrayMax(map.z)-1) return arr(map.z, arrayMax(map.z)-1, double) ;
  double xx = 0.01 * (x % 100) ;
  return (1-xx) * arr(map.z, xi, double) + xx * arr(map.z, xi+1, double) ;
}

/*****************************************/

typedef struct { 
  int lastPat[20], lastPos[20] ;
  double lastMap[20], glen[20] ;
  long nMinus[20], nPlus[20], len[20] ; 
} Hap4Stats ;
static Hap4Stats *stats ;

static double rateBoundary[20] = {
  0.1, 0.15, 0.2, 0.3, 0.5, 0.7,
  1.0, 1.5, 2.0, 3.0, 5.0, 7.0,
  10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 
  100.0, 1000.0
} ;

static void reportMinus (int varD, int x1, double g1, int x0, double g0)
{
  double rate = 1000000.0 * (g1 - g0) / (x1 - x0) ;
  int i = 0 ; while (rateBoundary[i] < rate) ++i ;
  ++stats[varD].nMinus[i] ;
  stats[varD].glen[i] += g1 - g0 ;
  stats[varD].len[i] += x1 - x0 ;
}

static void reportPlus (int varD, int x1, double g1, int x0, double g0)
{
  double rate = 1000000.0 * (g1 - g0) / (x1 - x0) ;
  int i = 0 ; while (rateBoundary[i] < rate) ++i ;
  ++stats[varD].nPlus[i] ;
  stats[varD].glen[i] += g1 - g0 ;
  stats[varD].len[i] += x1 - x0 ;
}

static void finalReport (void) 
{ 
  int v, i ;
  printf (" rate\tvar\t\tlen\tglen\tminus\t\tplus\n") ;
  for (v = 0 ; v < dictMax(variationDict) ; ++v)
    for (i = 0 ; i < 20 ; ++i)
      if (stats[v].nMinus[i] + stats[v].nPlus[i])
	printf ("%.2f\t%s\t%12ld\t%.4g\t%12ld\t%12ld\n", rateBoundary[i],
		dictName (variationDict, v), stats[v].len[i], stats[v].glen[i],
		stats[v].nMinus[i], stats[v].nPlus[i]) ;
}

void pbwt4hapsStats (PBWT *p)
{
  if (!p || !p->sites) die ("hap4stats called without a PBWT with sites") ;
  if (!map.x)
    { fprintf (logFile, "hap4stats called without a map - using a linear 1cM/Mb map\n") ;
      map.x = arrayCreate (2, int) ; 
      array(map.x,0,int) = arrp(p->sites,0,Site)->x ;
      array(map.x,1,int) = arrp(p->sites,arrayMax(p->sites)-1,Site)->x ;
      map.g = arrayCreate (2, double) ;
      array(map.g,0,double) = 0.0 ;
      array(map.g,1,int) = 0.000001 * (arr(map.x,1,int) - arr(map.x,0,int)) ;
      map.chrom = strdup (p->chrom) ;
      buildMap () ;
    }
  else if (strcmp (p->chrom, map.chrom)) 
    warn ("chrom mismatch in hap4stats: %s != %s", p->chrom, map.chrom) ;

  stats = mycalloc (dictMax(variationDict), Hap4Stats) ;

  int i, k ;
  int v ; 
  for (v = dictMax(variationDict) ; v-- ; ) 
    for (i = 0 ; i < p->M ; ++i) stats[v].lastPat[i] = -1 ;
  uchar *x = myalloc (p->M, uchar) ;
  
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  for (k = 0 ; k < p->N ; ++k)
    { int pos = arrp(p->sites, k, Site)->x ;
      double g = geneticMap (pos) ;
      v = arrp(p->sites,k,Site)->varD ;
      if (u->M - u->c >= 2)	/* at least two 1s */
	{ for (i = 0 ; i < p->M ; ++i) x[u->a[i]] = u->y[i] ;
	  for (i = 0 ; i+4 <= p->M ; i += 4)
	    if (x[i] + x[i+1] + x[i+2] + x[i+3] == 2) /* doubleton */
	      { int pat = x[i] + (x[i+1]<<1) + (x[i+2]<<2) + (x[i+3]<<3) ;
		if (stats[v].lastPat[i] >= 0)
		  { if (pat == stats[v].lastPat[i] || pat + stats[v].lastPat[i] == 15)
		      reportMinus (v, pos, g, 
				   stats[v].lastPos[i], stats[v].lastMap[i]) ;
		    else
		      reportPlus (v, pos, g, 
				  stats[v].lastPos[i], stats[v].lastMap[i]) ;
		  }
		stats[v].lastPat[i] = pat ; 
		stats[v].lastPos[i] = pos ; 
		stats[v].lastMap[i] = g ;
	      }
	}
      pbwtCursorForwardsRead (u) ;
    }

  finalReport () ;
  free (stats) ;
}

/************** end of file **************/

/*  File: pbwtMain.c
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
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr  2 08:05 2014 (rd)
 * Created: Thu Apr  4 12:05:20 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"

/*********************************************************/

static void prettyPlot (PBWT *p, FILE *fp, int K)
{
  int i, j ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  uchar **hap = pbwtHaplotypes (p) ;

  for (i = 0 ; i < K ; ++i)
    pbwtCursorForwardsRead (u) ;

  for (j = 0 ; j < p->M ; ++j)
    { for (i = K-100 ; i < K ; i++)
	putc (hap[u->a[j]][i]?'1':'0', fp) ;
      putc (' ',fp) ; putc (hap[u->a[j]][i++]?'1':'0', fp) ; putc (' ',fp) ;
      while (i < K+20)
	putc (hap[u->a[j]][i++]?'1':'0', fp) ;
      putc ('\n',fp) ;
    }
  pbwtCursorDestroy (u) ;
}

/*********************************************************/

static void exportSiteInfo (PBWT *p, FILE *fp, int f1, int f2)
/* print out d[] and y[] for sites with f1 <= f < f2 */
{
  int i, j, f, n = 0 ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;

  for (i = 0 ; i < p->N ; ++i)
    { f = p->M - u->c ;		/* number of 1s not 0s */
      if (f1 <= f && f < f2)	/* then print it out */
	{ for (j = 0 ; j < p->M ; ++j)
	    fprintf (fp, "%d %d ", u->y[j], i - u->d[j]) ;
	  fprintf (fp, "\n") ;
	  ++n ;
	}
      pbwtCursorForwardsReadAD (u, i) ;
    }
  pbwtCursorDestroy (u) ;
  fprintf (stderr, "%d rows exported with allele count f, %d <= f < %d\n", n, f1, f2) ;
}

/************ AF distribution **************/

static void siteFrequencySpectrum (PBWT *p)
{
  int i, j, n ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  Array hist = arrayCreate (p->M, int) ;
  int thresh[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9,
		   10, 20, 30, 40, 50, 60, 70, 80, 90,
		   100, 200, 300, 400, 500, 600, 700, 800, 900,
		   1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
		   10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
		   100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
		   1000000 } ;

  timeUpdate() ;

  for (i = 0 ; i < p->N ; ++i)
    { ++array(hist, p->M - u->c, int) ;
      pbwtCursorForwardsRead (u) ;
    }

  n = 0 ; j = 0 ;
  for (i = 1 ; i < p->M ; ++i)
    { n += array(hist,i,int) ;
      if (i == thresh[j])
	{ printf ("%d\t%d\n", thresh[j], n) ;
	  ++j ;
	  n = 0 ;
	}
    }
  printf ("%d\t%d\n", thresh[j], n) ;
}

/*********************************************************/

/* a couple of utilities for opening/closing files */

#define FOPEN(name,mode)  if (!strcmp (argv[1], "-")) fp = !strcmp(mode,"r") ? stdin : stdout ; else if (!(fp = fopen (argv[1],mode))) die ("failed to open %s file", name, argv[1])
#define FCLOSE if (strcmp(argv[1], "-")) fclose(fp)

int main (int argc, char *argv[])
{
  FILE *fp ;
  PBWT *p = 0 ;
  Array test ;

  pbwtInit () ;

  --argc ; ++argv ;
  if (!argc)			/* print help */
    { fprintf (stderr, "Usage: pbwt [ -<command> [options]* ]+\n") ;
      fprintf (stderr, "Commands:\n") ;
      fprintf (stderr, "  -check                    do various checks\n") ;
      fprintf (stderr, "  -stats                    print stats depending on commands; writes to stdout\n") ;
      fprintf (stderr, "  -read <file>              read pbwt file; '-' for stdin\n") ;
      fprintf (stderr, "  -readSites <file>         read sites file; '-' for stdin\n") ;
      fprintf (stderr, "  -readSamples <file>       read samples file; '-' for stdin\n") ;
      fprintf (stderr, "  -readMissing <file>       read missing file; '-' for stdin\n") ;
      fprintf (stderr, "  -readReverse <file>       read reverse file; '-' for stdin\n") ;
      fprintf (stderr, "  -readAll <rootname>       read .pbwt and if present .sites, .samples, .missing\n") ;
      fprintf (stderr, "  -readVcfGT <file>         read GTs from vcf or bcf file; '-' for stdin vcf only ; biallelic sites only - require diploid!\n") ;
      fprintf (stderr, "  -readVcfPL <file>         read PLs from vcf or bcf file; '-' for stdin vcf only ; biallelic sites only - require diploid!\n") ;
      fprintf (stderr, "  -readMacs <file>          read MaCS output file; '-' for stdin\n") ;
      fprintf (stderr, "  -readVcfq <file>          read VCFQ file; '-' for stdin\n") ;
      fprintf (stderr, "  -readGen <file> <chrom>   read impute2 gen file - must set chrom\n") ;
      fprintf (stderr, "  -readHap <file> <chrom>   read impute2 hap file - must set chrom\n") ;
      fprintf (stderr, "  -readPhase <file>         read Li and Stephens phase file\n") ;
      fprintf (stderr, "  -checkpoint <n>           checkpoint every n sites while reading\n") ;
      fprintf (stderr, "  -merge <file> ...         merge two or more pbwt files\n") ;
      fprintf (stderr, "  -write <file>             write pbwt file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeSites <file>        write sites file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeSamples <file>      write samples file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeMissing <file>      write missing file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeReverse <file>      write reverse file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeAll <rootname>      write .pbwt and if present .sites, .samples, .missing\n") ;
      fprintf (stderr, "  -writeImputeRef <rootname> write .imputeHaps and .imputeLegend\n") ;
      fprintf (stderr, "  -writeImputeHapsG <file>  write haplotype file for IMPUTE -known_haps_g\n") ;
      fprintf (stderr, "  -haps <file>              write haplotype file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeVcf <file>          write VCF file; '-' for stdout\n") ;
      fprintf (stderr, "  -subsites <fmin> <frac>   subsample <frac> sites with AF > <fmin>\n") ;
      fprintf (stderr, "  -subsample <start> <n>    subsample <n> samples from index <start>\n") ;
      fprintf (stderr, "  -subrange <start> <end>   cut down to sites in [start,end)\n") ;
      fprintf (stderr, "  -corruptSites <p> <q>     randomise fraction q of positions at fraction p of sites, according to site frequency\n") ;
      fprintf (stderr, "  -corruptSamples <p> <q>   randomise fraction q of positions for fraction p of samples, according to site frequency\n") ;
      fprintf (stderr, "  -selectSites <file>       select sites as in sites file\n") ;
      fprintf (stderr, "  -removeSites <file>       remove sites as in sites file\n") ;
      fprintf (stderr, "  -selectSamples <file>     select samples as in samples file\n") ;
      fprintf (stderr, "  -longWithin <L>           find matches within set longer than L\n") ;
      fprintf (stderr, "  -maxWithin                find maximal matches within set\n") ;
      fprintf (stderr, "  -matchNaive <file>        maximal match seqs in pbwt file to reference\n") ;
      fprintf (stderr, "  -matchIndexed <file>      maximal match seqs in pbwt file to reference\n") ;
      fprintf (stderr, "  -matchDynamic <file>      maximal match seqs in pbwt file to reference\n") ;
      fprintf (stderr, "  -imputeExplore <n>        n'th impute test\n") ;
      fprintf (stderr, "  -phase <k> <n>            phase with method k and n sparse pbwts\n") ;
      fprintf (stderr, "  -referencePhase <root>    phase current pbwt against reference whose root name is the argument - both pbwts need compatible sites! - only keep shared sites\n") ;
      fprintf (stderr, "  -referenceImpute <root>   impute current pbwt into reference whose root name is the argument - need compatible sites - does not rephase either pbwt\n") ;
      fprintf (stderr, "  -genotypeCompare <root>   compare genotypes with those from referencewhose root name is the argument - need compatible sites\n") ;
      fprintf (stderr, "  -imputeMissing            impute data marked as missing\n") ;
      fprintf (stderr, "  -paint                    output painting co-ancestry matrix\n") ;
      fprintf (stderr, "  -pretty <file> <k>        pretty plot at site k\n") ;
      fprintf (stderr, "  -sfs                      print site frequency spectrum (log scale)\n") ;
      fprintf (stderr, "  -siteInfo <file> <kmin> <kmax> export PBWT information at sites with allele count kmin <= k < kmax\n") ;
      fprintf (stderr, "  -buildReverse             build reverse pbwt\n") ;
    }

  timeUpdate() ;
  while (argc) {
    if (!(**argv == '-'))
      die ("not well formed command %s\nType pbwt without arguments for help", *argv) ;
    else if (!strcmp (argv[0], "-check"))
      { isCheck = TRUE ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-stats"))
      { isStats = TRUE ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-merge") && argc > 1)
    { 
        int i, nfiles = 0;
        const char **files = calloc(argc, sizeof(char*));
        for (i=1; i<argc; i++) 
        {
            if ( argv[i][0] == '-' ) break; // hopefully no one names their files to start with "-"
            files[nfiles++] = argv[i];
        }
        if ( nfiles>1 ) p = pbwtMerge(files, nfiles); 
        free(files);
        argc -= nfiles+1 ; argv += nfiles+1 ; 
    }
    else if (!strcmp (argv[0], "-haps") && argc > 1)
      { FOPEN("haps","w") ; pbwtWriteHaplotypes (fp, p) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-read") && argc > 1)
      { if (p) pbwtDestroy (p) ; FOPEN("read","r") ; p = pbwtRead (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readSites") && argc > 1)
      { FOPEN("readSites","r") ; pbwtReadSites (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readSamples") && argc > 1)
      { FOPEN("readSamples","r") ; pbwtReadSamples (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readMissing") && argc > 1)
      { FOPEN("readMissing","r") ; pbwtReadMissing (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readReverse") && argc > 1)
      { FOPEN("readReverse","r") ; pbwtReadReverse (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readAll") && argc > 1)
      { p = pbwtReadAll (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readVcfGT") && argc > 1)
      { if (p) pbwtDestroy (p) ; p = pbwtReadVcfGT (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readVcfPL") && argc > 1)
      { if (p) pbwtDestroy (p) ; p = pbwtReadVcfPL (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readMacs") && argc > 1)
      { if (p) pbwtDestroy (p) ; FOPEN("readMacs","r") ; p = pbwtReadMacs (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readVcfq") && argc > 1)
      { if (p) pbwtDestroy (p) ; FOPEN("readVcfq","r") ; p = pbwtReadVcfq (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readGen") && argc > 2)
      { if (p) pbwtDestroy (p) ; FOPEN("readGen","r") ; p = pbwtReadGen (fp, argv[2]) ; FCLOSE ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-readHap") && argc > 2)
    { if (p) pbwtDestroy (p) ; FOPEN("readHap","r") ; p = pbwtReadHap (fp, argv[2]) ; FCLOSE ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-readPhase") && argc > 1)
      { if (p) pbwtDestroy (p) ; FOPEN("readPhase","r") ; p = pbwtReadPhase (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-write") && argc > 1)
      { FOPEN("write","w") ; pbwtWrite (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeSites") && argc > 1)
      { FOPEN("writeSites","w") ; pbwtWriteSites (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeSamples") && argc > 1)
      { FOPEN("writeSamples","w") ; pbwtWriteSamples (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeMissing") && argc > 1)
      { FOPEN("writeMissing","w") ; pbwtWriteMissing (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeReverse") && argc > 1)
      { FOPEN("writeReverse","w") ; pbwtWriteReverse (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeAll") && argc > 1)
      { pbwtWriteAll (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeImputeRef") && argc > 1)
      { pbwtWriteImputeRef (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeImputeHapsG") && argc > 1)
      { FOPEN("writeImputeHaps","w") ; pbwtWriteImputeHapsG (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeGen") && argc > 1)
      { FOPEN("writeGen","w") ; pbwtWriteGen (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeVcf") && argc > 1)
      { pbwtWriteVcf (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-checkpoint") && argc > 1)
      { nCheckPoint = atoi (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-subsample") && argc > 2)
      { p = pbwtSubSampleInterval (p, atoi(argv[1]), atoi(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-selectSamples") && argc > 2)
      { FOPEN("selectSamples","r") ; p = pbwtSelectSamples (p, fp) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-subsites") && argc > 2)
      { p = pbwtSubSites (p, atof(argv[1]), atof(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-selectSites") && argc > 2)
      { FOPEN("selectSites","r") ; char *chr = 0 ; Array sites = pbwtReadSitesFile (fp, &chr) ;
	if (strcmp (chr, p->chrom)) die ("chromosome mismatch in selectSites") ;
	p = pbwtSelectSites (p, sites, FALSE) ; free (chr) ; arrayDestroy (sites) ;
	argc -= 2 ; argv += 2 ; 
      }
    else if (!strcmp (argv[0], "-removeSites") && argc > 2)
      { FOPEN("removeSites","r") ; char *chr = 0 ; Array sites = pbwtReadSitesFile (fp, &chr) ;
	if (p->chrom && strcmp (chr, p->chrom)) die ("chromosome mismatch in removeSites") ;
	p = pbwtRemoveSites (p, sites, FALSE) ; free (chr) ; arrayDestroy (sites) ;
	argc -= 2 ; argv += 2 ; 
      }
    else if (!strcmp (argv[0], "-subrange") && argc > 2)
      { p = pbwtSubRange (p, atoi(argv[1]), atoi(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-corruptSites") && argc > 2)
      { p = pbwtCorruptSites (p, atof(argv[1]), atof(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-corruptSamples") && argc > 2)
      { p = pbwtCorruptSamples (p, atof(argv[1]), atof(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-buildReverse"))
      { pbwtBuildReverse (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-pretty") && argc > 2)
      { FOPEN("prettyPlot","w") ; prettyPlot (p, fp, atoi(argv[2])) ; FCLOSE ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-siteInfo") && argc > 3)
      { FOPEN("siteInfo","w") ; exportSiteInfo (p, fp, atoi(argv[2]), atoi(argv[3])) ; FCLOSE ; argc -= 4 ; argv += 4 ; }
    else if (!strcmp (argv[0], "-sfs"))
      { siteFrequencySpectrum (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-maxWithin"))
      { pbwtLongMatches (p, 0) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-longWithin") && argc > 1)
      { pbwtLongMatches (p, atoi(argv[1])) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-matchNaive") && argc > 1)
      { FOPEN("matchNaive","r") ; matchSequencesNaive (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-matchIndexed") && argc > 1)
      { FOPEN("matchIndexed","r") ; matchSequencesIndexed (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-matchDynamic") && argc > 1)
      { FOPEN("matchDynamic","r") ; matchSequencesDynamic (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-imputeExplore") && argc > 1)
      { imputeExplore (p, atoi(argv[1])) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-phase") && argc > 2)
      { p = phase (p, atoi(argv[1]), atoi(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-referencePhase") && argc > 1)
      { p = referencePhase (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-referenceImpute") && argc > 1)
      { p = referenceImpute (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-genotypeCompare") && argc > 1)
      { genotypeCompare (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-imputeMissing"))
      { p = imputeMissing (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-paint"))
      { paintAncestryMatrix (p) ; argc -= 1 ; argv += 1 ; }
    else
      die ("unrecognised command %s\nType pbwt without arguments for help", *argv) ;
    timeUpdate() ;
  }
  if (p) pbwtDestroy(p) ;
  if (variationDict) dictDestroy(variationDict);
  sampleDestroy();
  fgetword (NULL) ;	// to keep valgrind happy, free malloced memory
  return 0 ;
}

/******************* end of file *******************/

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
 * Last edited: Dec 14 13:44 2015 (rd)
 * paintSparse added
 * Created: Thu Apr  4 12:05:20 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"
#include "version.h"

/*********************************************************/

#include "math.h"

static PBWT *playGround (PBWT *p)
{
  int k, j ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  double *d = 0, sumDiff2 = 0 ;

  for (k = 0 ; k < p->N ; ++k)
    { double psum = 0, xsum = 0, pxsum = 0 ;
      d = pbwtDosageRetrieve (p, u, d, k) ;
      for (j = 0 ; j < p->M ; ++j)
	{ psum += d[j] ;
	  if (d[j]) { xsum++ ; pxsum += d[j] ; }
	}
      psum /= p->M ; xsum /= p->M ; pxsum /= p->M ;
      double varProd = psum*(1.0-psum)*xsum*(1.0-xsum) ;
      double info = varProd ? (pxsum - psum*psum)/sqrt(varProd) : 1.0 ;
      double diff = (info - arrp(p->sites,k,Site)->imputeInfo) ;
      sumDiff2 += diff*diff ;
      pbwtCursorForwardsRead (u) ;
    }

  printf ("RMS info to zInfo %.4f\n", sqrt (sumDiff2/p->N)) ;

  pbwtCursorDestroy (u) ; free (d) ;
  return p ;
}

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
  fprintf (logFile, "%d rows exported with allele count f, %d <= f < %d\n", n, f1, f2) ;
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

  timeUpdate(logFile) ;

  FILE *fp ;
  if (p->sites) { fp = fopen ("sites.freq", "w") ; if (!fp) die ("can't open sites.freq") ; }

  for (i = 0 ; i < p->N ; ++i)
    { ++array(hist, p->M - u->c, int) ;
      if (p->sites)
	{ Site *s = arrp(p->sites,i,Site) ;
	  s->freq = 1.0 - (double)u->c/(double)p->M ;
	  fprintf (fp, "%s\t%d\t%.6f\t%s\n", p->chrom, s->x, s->freq, dictName (variationDict, s->varD)) ;
	}
      pbwtCursorForwardsRead (u) ;
    }
  if (p->sites) fclose (fp) ;

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

char *commandLine = "" ;

static void recordCommandLine (int argc, char *argv[])
 {
   if (!argc) return ;

   int i, len = 0 ;
   for (i = 0 ; i < argc ; ++i) len += (1 + strlen(argv[i])) ;
   commandLine = myalloc (len, char) ;
   strcpy (commandLine, argv[0]) ;
   for (i = 1 ; i < argc ; ++i) 
     { strcat (commandLine, " ") ; 
       strcat (commandLine, argv[i]) ;
     }
 }

/*********************************************************/

/* a couple of utilities for opening/closing files */

#define FOPEN(name,mode)  if (!strcmp (argv[1], "-")) fp = !strcmp(mode,"r") ? stdin : stdout ; else if (!(fp = fopen (argv[1],mode))) die ("failed to open %s file %s", name, argv[1])
#define FCLOSE if (strcmp(argv[1], "-")) fclose(fp)
#define LOPEN(name,mode)  if (!strcmp (argv[2], "-")) lp = !strcmp(mode,"r") ? stdin : stdout ; else if (!(lp = fopen (argv[2],mode))) die ("failed to open %s file", name, argv[2])
#define LCLOSE if (strcmp(argv[2], "-")) fclose(lp)
#define LOGOPEN(name) if (!strcmp (argv[1], "-")) logFile = stderr ; else if (!(logFile = fopen (argv[1],"w"))) die ("failed to open %s file %s", name, argv[1])
#define LOGCLOSE if (logFile && !(logFile==stderr)) fclose(logFile)

const char *pbwtCommitHash(void)
{
  return PBWT_COMMIT_HASH ;
}

FILE *logFile ; /* log file pointer */

int main (int argc, char *argv[])
{
  FILE *fp ;
  FILE *lp ;
  PBWT *p = 0 ;
  Array test ;
  char *referenceFasta = NULL;
  int isXY = 0 ;

  logFile = stderr ;

  pbwtInit () ;

  --argc ; ++argv ;
  recordCommandLine (argc, argv) ;

  if (!argc)			/* print help */
    { fprintf (stderr, "Program: pbwt\n") ;
      fprintf (stderr, "Version: %d.%d%s%s (using htslib %s)\n", pbwtMajorVersion, pbwtMinorVersion, 
                       strcmp(pbwtCommitHash(),"")==0 ? "" : "-", pbwtCommitHash(), pbwtHtslibVersionString()) ;
      fprintf (stderr, "Contact: Richard Durbin [rd@sanger.ac.uk]\n") ;
      fprintf (stderr, "Usage: pbwt [ -<command> [options]* ]+\n") ;
      fprintf (stderr, "Commands:\n") ;
      fprintf (stderr, "  -log <file>               log file; '-' for stderr\n") ;
      fprintf (stderr, "  -check                    do various checks\n") ;
      fprintf (stderr, "  -stats                    print stats depending on commands; writes to stdout\n") ;
      fprintf (stderr, "  -read <file>              read pbwt file; '-' for stdin\n") ;
      fprintf (stderr, "  -readSites <file>         read sites file; '-' for stdin\n") ;
      fprintf (stderr, "  -readSamples <file>       read samples file; '-' for stdin\n") ;
      fprintf (stderr, "  -readMissing <file>       read missing file; '-' for stdin\n") ;
      fprintf (stderr, "  -readDosage <file>        read dosage file; '-' for stdin\n") ;
      fprintf (stderr, "  -readReverse <file>       read reverse file; '-' for stdin\n") ;
      fprintf (stderr, "  -readAll <rootname>       read .pbwt and if present .sites, .samples, .missing - note not by default dosage\n") ;
      fprintf (stderr, "  -loadSamples <file>       load sample metadata from a FMF, IMPUTE2 or FAM samples file; '-' for stdin\n") ;
      fprintf (stderr, "  -readVcfGT <file>         read GTs from vcf or bcf file; '-' for stdin vcf only\n") ;
      fprintf (stderr, "  -readVcfPL <file>         read PLs from vcf or bcf file; '-' for stdin vcf only\n") ;
      fprintf (stderr, "  -readMacs <file>          read MaCS output file; '-' for stdin\n") ;
      fprintf (stderr, "  -readVcfq <file>          read VCFQ file; '-' for stdin\n") ;
      fprintf (stderr, "  -readGen <file> <chrom>   read impute2 gen file - must set chrom\n") ;
      fprintf (stderr, "  -readHap <file> <chrom>   read impute2 hap file - must set chrom\n") ;
      fprintf (stderr, "  -readHapLegend <hap_file> <legend_file> <chrom>\n") ;
      fprintf (stderr, "                            read impute2 hap and legend file - must set chrom\n") ;
      fprintf (stderr, "  -readPhase <file>         read Li and Stephens phase file\n") ;
      fprintf (stderr, "  -checkpoint <n>           checkpoint every n sites while reading\n") ;
      fprintf (stderr, "  -merge <file> ...         merge two or more pbwt files\n") ;
      fprintf (stderr, "  -write <file>             write pbwt file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeSites <file>        write sites file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeSamples <file>      write samples file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeMissing <file>      write missing file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeDosage <file>       write missing file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeReverse <file>      write reverse file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeAll <rootname>      write .pbwt and if present .sites, .samples, .missing, .dosage\n") ;
      fprintf (stderr, "  -writeImputeRef <rootname> write .imputeHaps and .imputeLegend\n") ;
      fprintf (stderr, "  -writeImputeHapsG <file>  write haplotype file for IMPUTE -known_haps_g\n") ;
      fprintf (stderr, "  -writePhase <file>        write FineSTRUCTURE/ChromoPainter input format (Impute/ShapeIT output format) phase file\n") ;
      fprintf (stderr, "  -writeTransposeHaplotypes <file>    write transposed haplotype file (one hap per row); '-' for stdout\n") ;
      fprintf (stderr, "  -haps <file>              write haplotype file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeGen <file>          write impute2 gen file; '-' for stdout\n") ;
      fprintf (stderr, "  -writeVcf|-writeVcfGz|-writeBcf|-writeBcfGz <file>\n") ;
      fprintf (stderr, "                            write VCF or BCF; uncompressed or bgzip (Gz) compressed file; '-' for stdout\n") ;
      fprintf (stderr, "  -referenceFasta <file>    reference fasta filename for VCF/BCF writing (optional)\n") ;
      fprintf (stderr, "  -X                        treat males as haploid, females diploid as on non-PAR chrX\n") ;
      fprintf (stderr, "  -Y                        treat males as haploid, skip females as on chrY\n") ;
      fprintf (stderr, "  -subsites <fmin> <frac>   subsample <frac> sites with AF > <fmin>\n") ;
      fprintf (stderr, "  -subsample <start> <n>    subsample <n> samples from index <start>\n") ;
      fprintf (stderr, "  -subrange <start> <end>   cut down to sites in [start,end)\n") ;
      fprintf (stderr, "  -corruptSites <p> <q>     randomise fraction q of positions at fraction p of sites, according to site frequency\n") ;
      fprintf (stderr, "  -corruptSamples <p> <q>   randomise fraction q of positions for fraction p of samples, according to site frequency\n") ;
      fprintf (stderr, "  -copySamples <M> <len>    make M new samples copied from current haplotypes with mean switch length len\n") ;
      fprintf (stderr, "  -selectSites <file>       select sites as in sites file\n") ;
      fprintf (stderr, "  -removeSites <file>       remove sites as in sites file\n") ;
      fprintf (stderr, "  -selectSamples <file>     select samples as in samples file\n") ;
      fprintf (stderr, "  -removeSamples <file>     remove samples as in samples file\n") ;
      fprintf (stderr, "  -longWithin <L>           find matches within set longer than L\n") ;
      fprintf (stderr, "  -maxWithin                find maximal matches within set\n") ;
      fprintf (stderr, "  -matchNaive <file>        maximal match seqs in pbwt file to reference\n") ;
      fprintf (stderr, "  -matchIndexed <file>      maximal match seqs in pbwt file to reference\n") ;
      fprintf (stderr, "  -matchDynamic <file>      maximal match seqs in pbwt file to reference\n") ;
      fprintf (stderr, "  -imputeExplore <n>        n'th impute test\n") ;
      fprintf (stderr, "  -phase <n>                phase with n sparse pbwts\n") ;
      fprintf (stderr, "  -referencePhase <root>    phase current pbwt against reference whose root name is the argument - only keeps shared sites\n") ;
      fprintf (stderr, "  -referenceImpute <root> [nSparse=1] [fSparse=1]  impute current pbwt into reference whose root name is the first argument;\n") ;
      fprintf (stderr, "                            does not rephase either pbwt; optional nSparse > 1 also does sparse matching, fSparse is relative weight\n") ;
      fprintf (stderr, "  -genotypeCompare <root>   compare genotypes with those from reference whose root name is the argument - need compatible sites\n") ;
      fprintf (stderr, "  -imputeMissing            impute data marked as missing\n") ;
      fprintf (stderr, "  -fitAlphaBeta <model>     fit probabilistic model 1..3\n") ;
      fprintf (stderr, "  -llCopyModel <theta> <rho>  log likelihood of Li-Stephens model\n") ;
      fprintf (stderr, "  -paint <fileNameRoot> [n] output painting co-ancestry matrix to fileroot, optionally specififying the number per region\n") ;
      fprintf (stderr, "  -paintSparse <fileNameRoot> [n] output sparse painting to fileroot, optionally specififying the number per region\n") ;
      fprintf (stderr, "  -pretty <file> <k>        pretty plot at site k\n") ;
      fprintf (stderr, "  -sfs                      print site frequency spectrum (log scale) - also writes sites.freq file\n") ;
      fprintf (stderr, "  -refFreq <file>           read site frequency information into the refFreq field of current sites\n") ;
      fprintf (stderr, "  -siteInfo <file> <kmin> <kmax> export PBWT information at sites with allele count kmin <= k < kmax\n") ;
      fprintf (stderr, "  -buildReverse             build reverse pbwt\n") ;
      fprintf (stderr, "  -readGeneticMap <file>    read Oxford format genetic map file\n") ;
      fprintf (stderr, "  -4hapsStats               mu:rho 4 hap test stats\n") ;
    }

  timeUpdate(logFile) ;
  while (argc) {
    if (!(**argv == '-'))
      die ("not well formed command %s\nType pbwt without arguments for help", *argv) ;
    else if (!strcmp (argv[0], "-check"))
      { isCheck = TRUE ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-stats"))
      { isStats = TRUE ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-X"))
      { isXY = 2 ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-Y"))
      { isXY = 1 ; argc -= 1 ; argv += 1 ; }
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
    else if (!strcmp (argv[0], "-log") && argc > 1)
      { LOGCLOSE ; LOGOPEN("log") ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-haps") && argc > 1)
      { FOPEN("haps","w") ; pbwtWriteHaplotypes (fp, p) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-read") && argc > 1)
      { if (p) pbwtDestroy (p) ; FOPEN("read","r") ; p = pbwtRead (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readSites") && argc > 1)
      { FOPEN("readSites","r") ; pbwtReadSites (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readSamples") && argc > 1)
      { FOPEN("readSamples","r") ; pbwtReadSamples (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-loadSamples") && argc > 1)
      { FOPEN("loadSamples","r") ; pbwtReadSamplesFile (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readMissing") && argc > 1)
      { FOPEN("readMissing","r") ; pbwtReadMissing (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readDosage") && argc > 1)
      { FOPEN("readMissing","r") ; pbwtReadDosage (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readReverse") && argc > 1)
      { FOPEN("readReverse","r") ; pbwtReadReverse (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readAll") && argc > 1)
      { p = pbwtReadAll (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-readVcfGT") && argc > 1)
      { if (p) pbwtDestroy (p) ; p = pbwtReadVcfGT (argv[1], isXY) ; argc -= 2 ; argv += 2 ; }
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
    else if (!strcmp (argv[0], "-readHapLegend") && argc > 3)
    { if (p) pbwtDestroy (p) ; FOPEN("readHap","r") ; LOPEN("readHap","r") ; p = pbwtReadHapLegend (fp, lp, argv[3]) ; FCLOSE ; LCLOSE ; argc -= 4 ; argv += 4 ; }
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
    else if (!strcmp (argv[0], "-writeDosage") && argc > 1)
      { FOPEN("writeMissing","w") ; pbwtWriteDosage (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
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
    else if (!strcmp (argv[0], "-writePhase") && argc > 1)
      { pbwtWritePhase (p,argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeTransposedHaplotypes") && argc > 1)
      { FOPEN("writeTransposedHaplotypes",argv[0]) ; pbwtWriteTransposedHaplotypes (p, fp) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-referenceFasta") && argc > 1)
      { referenceFasta = strdup(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeVcf") && argc > 1)
      { pbwtWriteVcf (p, argv[1], referenceFasta, "w") ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeVcfGz") && argc > 1)
      { pbwtWriteVcf (p, argv[1], referenceFasta, "wz") ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeBcf") && argc > 1)
      { pbwtWriteVcf (p, argv[1], referenceFasta, "wbu") ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-writeBcfGz") && argc > 1)
      { pbwtWriteVcf (p, argv[1], referenceFasta, "wb") ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-checkpoint") && argc > 1)
      { nCheckPoint = atoi (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-subsample") && argc > 2)
      { p = pbwtSubSampleInterval (p, atoi(argv[1]), atoi(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-selectSamples") && argc > 2)
      { FOPEN("selectSamples","r") ; p = pbwtSelectSamples (p, fp) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-removeSamples") && argc > 2)
      { FOPEN("removeSamples","r") ; p = pbwtRemoveSamples (p, fp) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-subsites") && argc > 2)
      { p = pbwtSubSites (p, atof(argv[1]), atof(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-selectSites") && argc > 1)
      { FOPEN("selectSites","r") ; char *chr = 0 ; Array sites = pbwtReadSitesFile (fp, &chr) ;
	if (strcmp (chr, p->chrom)) die ("chromosome mismatch in selectSites") ;
	p = pbwtSelectSites (p, sites, FALSE) ; free (chr) ; arrayDestroy (sites) ;
	argc -= 2 ; argv += 2 ; 
      }
    else if (!strcmp (argv[0], "-removeSites") && argc > 1)
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
    else if (!strcmp (argv[0], "-copySamples") && argc > 2)
      { p = pbwtCopySamples (p, atoi(argv[1]), atof(argv[2])) ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-buildReverse"))
      { pbwtBuildReverse (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-pretty") && argc > 2)
      { FOPEN("prettyPlot","w") ; prettyPlot (p, fp, atoi(argv[2])) ; FCLOSE ; argc -= 3 ; argv += 3 ; }
    else if (!strcmp (argv[0], "-siteInfo") && argc > 3)
      { FOPEN("siteInfo","w") ; exportSiteInfo (p, fp, atoi(argv[2]), atoi(argv[3])) ; FCLOSE ; argc -= 4 ; argv += 4 ; }
    else if (!strcmp (argv[0], "-sfs"))
      { siteFrequencySpectrum (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-refFreq") && argc > 1)
      { FOPEN("refFreq","r") ; pbwtReadRefFreq (p, fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
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
    else if (!strcmp (argv[0], "-phase") && argc > 1)
      { p = phase (p, atoi(argv[1])) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-referencePhase") && argc > 1)
      { p = referencePhase (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-referenceImpute") && argc > 1)
      { int nSparse = 1 ; double fSparse = 1.0 ;
	char *fileNameRoot = argv[1] ;  argc -= 2 ; argv += 2 ; 
	if (argc && argv[0][0] != '-')
	  { if (!(nSparse = atoi(argv[0]))) die ("bad refImpute nSparse %s", argv[0]) ;
	    else { --argc ; ++argv ; }
	  }
	if (argc && argv[0][0] != '-')
	  { if (!(fSparse = atof(argv[0]))) die ("bad refImpute fSparse %s", argv[0]) ;
	    else { --argc ; ++argv ; }
	  }
	p = referenceImpute (p, fileNameRoot, nSparse, fSparse) ;
      }
    else if (!strcmp (argv[0], "-genotypeCompare") && argc > 1)
      { genotypeCompare (p, argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-imputeMissing"))
      { p = imputeMissing (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-fitAlphaBeta") && argc > 1)
      { pbwtFitAlphaBeta (p, atoi(argv[1])) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-llCopyModel") && argc > 2)
      { pbwtLogLikelihoodCopyModel (p, atof(argv[1]), atof(argv[2])) ; 
	argc -= 3 ; argv += 3 ; 
      }
    else if (!strcmp (argv[0], "-readGeneticMap") && argc > 1)
      {  FOPEN("readGeneticMap","r") ; readGeneticMap (fp) ; FCLOSE ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (argv[0], "-4hapsStats"))
      { pbwt4hapsStats (p) ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (argv[0], "-paint") && argc > 1)
      { 
	int npr=100;
       	if(argc>2) if(argv[2][0] !='-') npr=atoi(argv[2]);
	paintAncestryMatrix (p, argv[1],npr) ; argc -= 2 ; argv += 2 ; 
       	if(argc>0) if(argv[0][0] !='-') {--argc;++argv; }
      }
    else if (!strcmp (argv[0], "-paintSparse") && argc > 1)
      { 
	int npr=100;
       	if(argc>2) if(argv[2][0] !='-') npr=atoi(argv[2]);
	paintAncestryMatrixSparse (p, argv[1],npr,0) ; argc -= 2 ; argv += 2 ; 
       	if(argc>0) if(argv[0][0] !='-') {--argc;++argv; }
      }
    else if (!strcmp (argv[0], "-play"))
      { p = playGround (p) ; argc -= 1 ; argv += 1 ; }
    else
      die ("unrecognised command %s\nType pbwt without arguments for help", *argv) ;
    timeUpdate(logFile) ;
  }
  if (p) pbwtDestroy(p) ;
  if (variationDict) dictDestroy(variationDict);
  sampleDestroy();
  metaDataDestroy();
  if (*commandLine) free(commandLine);
  fgetword (NULL) ;	// to keep valgrind happy, free malloced memory
  LOGCLOSE ;
  return 0 ;
}

/******************* end of file *******************/

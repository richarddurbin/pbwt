/*  File: pbwtHtslib.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013
 *-------------------------------------------------------------------
 * Description: all the pbwt stuff that uses htslib, e.g. reading/writing vcf or bcf files
 * Exported functions:
 * HISTORY:
 * Last edited: Oct  3 15:48 2014 (rd)
 * * Sep 22 23:03 2014 (rd): change for 64bit arrays
 * Created: Thu Oct 17 12:20:04 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "pbwt.h"
#include <htslib/synced_bcf_reader.h>

static void readVcfSamples (PBWT *p, bcf_hdr_t *hr)
{
  int i, k ;

  p->samples = arrayCreate (p->M, int) ;
  for (i = 0 ; i < p->M/2 ; ++i)
    { int k = sampleAdd (hr->samples[i],0,0,0) ;
      array(p->samples, 2*i, int) = k ; /* assume diploid - could be cleverer */
      array(p->samples, 2*i+1, int) = k ;
    }
}

static int variation (PBWT *p, const char *ref, const char *alt)
{
  static char *buf = 0 ;
  static int buflen = 0 ;
  if (!buf) { buflen = 64 ; buf = myalloc (buflen, char) ; }
  int var ;
  if (strlen (ref) + strlen (alt) + 2 > buflen) 
    { do buflen *= 2 ; while (strlen (ref) + strlen (alt) + 2 > buflen) ;
      free (buf) ; buf = myalloc (buflen, char) ;
    }
  sprintf (buf, "%s\t%s", ref, alt) ;
  dictAdd (variationDict, buf, &var) ;
  return var ;
}

PBWT *pbwtReadVcfGT (char *filename)  /* read GTs from vcf/bcf using htslib */
{
  int i, j ;

  bcf_srs_t *sr = bcf_sr_init() ;
  if (!bcf_sr_add_reader (sr, filename)) die ("failed to open good vcf file\n") ;

  bcf_hdr_t *hr = sr->readers[0].header ;
  PBWT *p = pbwtCreate (bcf_hdr_nsamples(hr)*2, 0) ; /* assume diploid! */
  readVcfSamples (p, hr) ;
  p->sites = arrayCreate (10000, Site) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  uchar *x = myalloc (p->M, uchar) ;

  uchar *xMissing = myalloc (p->M+1, uchar) ;
  xMissing[p->M] = Y_SENTINEL ;  /* needed for efficient packing */
  long nMissing = 0 ;
  int nMissingSites = 0 ; 

  int mgt_arr = 0, *gt_arr = NULL ;
  while (bcf_sr_next_line (sr)) 
    { bcf1_t *line = bcf_sr_get_line(sr,0) ;
      const char* chrom = bcf_seqname(hr,line) ;
      if (!p->chrom) p->chrom = strdup (chrom) ;
      else if (strcmp (chrom, p->chrom)) break ;
      int pos = line->pos + 1 ;       // bcf coordinates are 0-based
      // if (line->n_allele != 2) continue ;  // not a biallelic site - skip these
      const char *ref = line->d.allele[0] ;

      // get a copy of GTs
      int ngt = bcf_get_genotypes(hr, line, &gt_arr, &mgt_arr) ;
      if (ngt <= 0) continue ;  // it seems that -1 is used if GT is not in the FORMAT
      if (ngt != p->M) die ("%d != %d GT values at %s:%d - not diploid?", 
          ngt, p->M, chrom, pos) ;

      memset (xMissing, 0, p->M) ;
      long wasMissing = nMissing ;

      /* copy the genotypes into array x[] */
      for (i = 0 ; i < p->M ; i++)
        { if (gt_arr[i] == bcf_int32_vector_end) 
            x[i] = bcf_gt_allele(gt_arr[i-1]); // treat haploid genotypes as diploid homozygous A/A
          if (gt_arr[i] == bcf_gt_missing)
            { x[i] = 0 ; /* use ref for now */
              xMissing[i] = 1 ;
              ++nMissing ;
            }
          else 
            x[i] = bcf_gt_allele(gt_arr[i]) ;  // convert from BCF binary to 0 or 1
        }

      /* split into biallelic sites filling in as REF ALT alleles */
      /* not in the REF/ALT site */
      for (i = 1 ; i < line->n_allele ; i++)
        {
          const char *alt = line->d.allele[i] ;

          /* and pack them into the PBWT */
          for (j = 0 ; j < p->M ; ++j) u->y[j] = x[u->a[j]] == i ? 1 : 0;
          pbwtCursorWriteForwards (u) ;

          /* store missing information, if there was any */
          if (nMissing > wasMissing)
            { if (!wasMissing)
                { p->zMissing = arrayCreate (10000, uchar) ;
                  array(p->zMissing, 0, uchar) = 0 ; /* needed so missing[] has offset > 0 */
                  p->missingOffset = arrayCreate (1024, long) ;
                }
              array(p->missingOffset, p->N, long) = arrayMax(p->zMissing) ;
              pack3arrayAdd (xMissing, p->M, p->zMissing) ; /* NB original order, not pbwt sort */
              nMissingSites++ ;
            }
          else if (nMissing)
            array(p->missingOffset, p->N, int) = 0 ;

          // add the site
          Site *s = arrayp(p->sites, p->N++, Site) ;
          s->x = pos ;
          s->varD = variation (p, ref, alt) ;          
        }

      if (nCheckPoint && !(p->N % nCheckPoint))  pbwtCheckPoint (u, p) ;
    }
  pbwtCursorToAFend (u, p) ;

  if (gt_arr) free (gt_arr) ;
  bcf_sr_destroy (sr) ;
  free (x) ; pbwtCursorDestroy (u) ;  
  free (xMissing) ;

  fprintf (stderr, "read genotypes from %s\n", filename) ;
  if (p->missingOffset) fprintf (stderr, "%ld missing values at %d sites\n", 
         nMissing, nMissingSites) ;

  return p ;
}

PBWT *pbwtReadVcfPL (char *filename)  /* read PLs from vcf/bcf using htslib */
{
  PBWT *p ;
  int i, j, k = 0 ;

  bcf_srs_t *sr = bcf_sr_init() ;
  if (!bcf_sr_add_reader (sr, filename)) die ("failed to open good vcf file\n") ;

  bcf_hdr_t *hr = sr->readers[0].header ;
  p = pbwtCreate (bcf_hdr_nsamples(hr)*2, 0) ; /* assume diploid! */
  readVcfSamples (p, hr) ;

  int mpl_arr = 0, *pl_arr = NULL;
  while (bcf_sr_next_line (sr)) 
    { ++k ;
      bcf1_t *line = bcf_sr_get_line(sr,0) ;
      const char* chrom = bcf_seqname(hr,line) ;
      int pos = line->pos + 1 ;       // bcf coordinates are 0-based
      if (line->n_allele != 2) continue ;  // not a biallelic site
      const char *ref = line->d.allele[0] ;
      const char *alt = line->d.allele[1] ;

      if (k <= 10) printf ("%s:%d %s %s", chrom, pos, ref, alt) ;

      // get a copy of the PL vectors
      int npl = bcf_get_format_int32(hr, line, "PL", &pl_arr, &mpl_arr) ;
      if (npl)
        { npl /= bcf_hdr_nsamples(hr) ;  // number of values per samples
          if (npl != 3) die ("%s:%d not a diploid site", chrom, pos) ; // not diploid
          for (i = 0 ; i < bcf_hdr_nsamples(hr) ; i++) // iterate over samples
            {
              for (j = 0 ; j < npl ; j++) // iterate over PL values (genotypes)
                {
                  // check for shorter vectors (haploid genotypes amongst diploid genotypes)
                  if (pl_arr[i*npl+j] == bcf_int32_vector_end) break ;
                  // skip missing values
                  if (pl_arr[i*npl+j] == bcf_int32_missing) continue ;

                  // do something with j-th PL value
                  if (k <= 10 && i < 10) printf ("%c%d", j?'.':' ', pl_arr[i*npl+j]) ;
                }
            }
        }
      if (k <= 10) printf("\n");
    }

  if (pl_arr) free (pl_arr) ;
  bcf_sr_destroy (sr) ;
  
  return p ;
}

void pbwtWriteVcf (PBWT *p, char *filename)  /* write vcf/bcf using htslib */
{
  htsFile *bcf_fp = NULL ;
  bcf_hdr_t *bcf_hdr = NULL ;

  // wb .. compressed BCF
  // wbu .. uncompressed BCF
  // wz .. compressed VCF
  // w .. uncompressed VCF
  // write out uncompressed VCF for the moment. BCF needs ##contig metadata complete in header
  bcf_fp = hts_open(filename,"w");
  if (!p) die ("pbwtWriteVcf called without a valid pbwt") ;
  if (!p->sites) die ("pbwtWriteVcf called without sites") ;
  int samples_known = 1 ;
  if (!p->samples)
    {
      fprintf (stderr, "Warning: pbwtWriteVcf called without samples... using fake sample names PBWT0, PBWT1 etc...\n") ;
      samples_known = 0 ;
    }

  // write header
  bcf_hdr = bcf_hdr_init("w") ;
  kstring_t str = {0,0,0} ;
  ksprintf(&str, "##pbwtVersion=%d.%d+htslib-%s", 
	   pbwtMajorVersion, pbwtMinorVersion, hts_version()) ;
  bcf_hdr_append(bcf_hdr, str.s) ;
  free(str.s) ;
  bcf_hdr_append(bcf_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">") ;
  bcf_hdr_append(bcf_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">") ;
  bcf_hdr_append(bcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">") ;
  
  int i, j ;
  for (i = 0 ; i < p->M/2 ; ++i)
    {
      if (samples_known)
        bcf_hdr_add_sample(bcf_hdr, sampleName(sample (p, 2*i))) ;
      else
        {
          kstring_t sname = {0,0,0} ;
          ksprintf(&sname, "PBWT%d", i) ;
          bcf_hdr_add_sample(bcf_hdr, sname.s) ;
          free(sname.s) ;
        }
    }
  bcf_hdr_add_sample(bcf_hdr, 0) ; /* required to update internal structures */
  bcf_hdr_write(bcf_fp, bcf_hdr) ;

  bcf1_t *bcf_rec = bcf_init1() ;
  uchar *hap = myalloc (p->M, uchar) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;

  for (i = 0 ; i < p->N ; ++i)
    {
      kstring_t rec = {0,0,0} ;
      Site *s = arrp(p->sites, i, Site) ;
      ksprintf(&rec, "%s\t%d\t.\t%s\t.\tPASS\t.\tGT", 
        p->chrom ? p->chrom : ".", 
        s->x, dictName(variationDict, s->varD));

      for (j = 0 ; j < p->M ; ++j)
        {
          hap[u->a[j]] = u->y[j] ;
        }
      int ac[2] = {0,0};
      for (j = 0 ; j < p->M ; j+=2)
        {
          // write out phased or unphased? phased flag could be stored in pbwt struct?
          // could add PS or other FORMAT fields here?
          // todo: handle missing data
          ksprintf(&rec, "\t%d|%d", hap[j], hap[j+1]) ;
          ac[hap[j]]++ ;
          ac[hap[j+1]]++ ;
        }
      int an = ac[0] + ac[1] ;

      // parse VCF record string into bcf record. 
      // could be quicker to load directly into bcf, but this is easy

      vcf_parse(&rec, bcf_hdr, bcf_rec) ; 

      // example of adding INFO fields
      bcf_update_info_int32(bcf_hdr, bcf_rec, "AC", &ac[1], 1) ;
      bcf_update_info_int32(bcf_hdr, bcf_rec, "AN", &an, 1) ;

      //write and progress
      bcf_write1(bcf_fp, bcf_hdr, bcf_rec) ;
      bcf_clear1(bcf_rec) ;

      free(rec.s) ;
      pbwtCursorForwardsRead(u) ;
    }

  // cleanup
  free(hap) ;
  /*  pbwtCursorDestroy(u) ; */
  bcf_hdr_destroy(bcf_hdr) ;
  bcf_destroy1(bcf_rec);
  hts_close(bcf_fp) ;
}

/******* end of file ********/

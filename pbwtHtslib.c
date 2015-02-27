/*  File: pbwtHtslib.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013
 *-------------------------------------------------------------------
 * Description: all the pbwt stuff that uses htslib, e.g. reading/writing vcf or bcf files
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 13 15:17 2014 (rd)
 * * Sep 22 23:03 2014 (rd): change for 64bit arrays
 * Created: Thu Oct 17 12:20:04 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "pbwt.h"
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>

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
      char *ref, *REF; 
      ref = REF = strdup(line->d.allele[0]);
      while ( *ref++ = toupper(*ref) );

      // get a copy of GTs
      int ngt = bcf_get_genotypes(hr, line, &gt_arr, &mgt_arr) ;
      if (ngt <= 0) continue ;  // it seems that -1 is used if GT is not in the FORMAT
      if (ngt != p->M && p->M != 2*ngt) die ("%d != %d GT values at %s:%d - not haploid or diploid?", 
          ngt, p->M, chrom, pos) ;

      memset (xMissing, 0, p->M) ;
      long wasMissing = nMissing ;
      /* copy the genotypes into array x[] */
      if (p->M == 2*ngt) // all GTs haploid: treat haploid genotypes as diploid homozygous A/A
        {
          for (i = 0 ; i < ngt ; i++)
            { if (gt_arr[i] == bcf_gt_missing)
                { x[2*i] = 0 ;
                  x[2*i+1] = 0; /* use ref for now */
                  xMissing[2*i] = 1 ;
                  xMissing[2*i+1] = 1;
                  nMissing+=2 ;
                }
              else {
                x[2*i] = bcf_gt_allele(gt_arr[i]) ;  // convert from BCF binary to 0 or 1
                x[2*i+1] = x[2*i] ;  // convert from BCF binary to 0 or 1
              }
            }
        }
      else
        {
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
        }

      BOOL no_alt = line->n_allele == 1;
      int n_allele = no_alt ? 2 : line->n_allele;

      /* split into biallelic sites filling in as REF ALT alleles */
      /* not in the REF/ALT site */
      for (i = 1 ; i < n_allele ; i++)
        {
          char *alt, *ALT; 
          alt = ALT = no_alt ? "." : strdup(line->d.allele[i]) ;
          if (!no_alt) while ( *alt++ = toupper(*alt) );

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
            array(p->missingOffset, p->N, long) = 0 ;

          // add the site
          Site *s = arrayp(p->sites, p->N++, Site) ;
          s->x = pos ;
          s->varD = variation (p, REF, ALT) ;          
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

static void pbwtSetContigs(bcf_hdr_t *hdr, faidx_t *fai)
{
  int i, n = faidx_nseq(fai) ;
  for (i=0; i<n; i++)
    {
      const char *seq = faidx_iseq(fai,i) ;
      int len = faidx_seq_len(fai, seq) ;
      bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq, len) ;
    }
}

void pbwtWriteVcf (PBWT *p, char *filename, char *referenceFasta, char *mode)
{
  htsFile *fp = NULL ;
  bcf_hdr_t *bcfHeader = NULL ;

  fp = hts_open(filename,mode) ;
  if (!fp) die ("could not open file for writing: %s", filename) ;
  if (!p) die ("pbwtWriteVcf called without a valid pbwt") ;
  if (!p->sites) die ("pbwtWriteVcf called without sites") ;
  if (!p->samples) fprintf (stderr, "Warning: pbwtWriteVcf called without samples... using fake sample names PBWT0, PBWT1 etc...\n") ;

  // write header
  bcfHeader = bcf_hdr_init("w") ;
  if (referenceFasta)
    {
      faidx_t *faidx = fai_load(referenceFasta);
      if ( !faidx ) die ("Could not load the reference %s. Has the fasta been indexed with 'samtools faidx'?\n", referenceFasta);
      pbwtSetContigs(bcfHeader, faidx);
      fai_destroy(faidx);
    }
  else if (p->chrom)
    {
      bcf_hdr_printf(bcfHeader, "##contig=<ID=%s,length=%d>", p->chrom, 0x7fffffff);   // MAX_CSI_COOR
    }
  kstring_t str = {0,0,0} ;
  ksprintf(&str, "##pbwtVersion=%d.%d+htslib-%s", 
	   pbwtMajorVersion, pbwtMinorVersion, hts_version()) ;
  bcf_hdr_append(bcfHeader, str.s) ;
  free(str.s) ;
  bcf_hdr_append(bcfHeader, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">") ;
  bcf_hdr_append(bcfHeader, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">") ;
  bcf_hdr_append(bcfHeader, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">") ;
  
  int i, j ;
  for (i = 0 ; i < p->M/2 ; ++i)
    {
      if (p->samples)
        bcf_hdr_add_sample(bcfHeader, sampleName(sample (p, 2*i))) ;
      else
        {
          kstring_t sname = {0,0,0} ;
          ksprintf(&sname, "PBWT%d", i) ;
          bcf_hdr_add_sample(bcfHeader, sname.s) ;
          free(sname.s) ;
        }
    }
  bcf_hdr_add_sample(bcfHeader, 0) ; /* required to update internal structures */
  bcf_hdr_write(fp, bcfHeader) ;

  bcf1_t *bcfRecord = bcf_init1() ;
  uchar *hap = myalloc (p->M, uchar) ;
  int32_t *gts = myalloc (p->M, int32_t) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;

  for (i = 0 ; i < p->N ; ++i)
    {
      Site *s = arrp(p->sites, i, Site) ;
      bcf_float_set_missing(bcfRecord->qual) ;
      bcfRecord->rid = bcf_hdr_name2id(bcfHeader, p->chrom) ;
      bcfRecord->pos = s->x - 1 ;
      char *als = strdup( dictName(variationDict, s->varD) ), *ss = als ;
      while ( *ss ) { if ( *ss=='\t' ) *ss = ',' ; ss++ ; }
      bcf_update_alleles_str(bcfHeader, bcfRecord, als) ;
      bcf_add_filter(bcfHeader, bcfRecord, bcf_hdr_id2int(bcfHeader, BCF_DT_ID, "PASS")) ;

      for (j = 0 ; j < p->M ; ++j)
        {
          hap[u->a[j]] = u->y[j] ;
        }
      int ac[2] = {0,0};
      for (j = 0 ; j < p->M ; j+=2)
        {
          // todo: handle missing data
          gts[j] = bcf_gt_unphased(hap[j]);
          gts[j+1] = bcf_gt_phased(hap[j+1]);
          ac[hap[j]]++ ;
          ac[hap[j+1]]++ ;
        }
      int an = ac[0] + ac[1] ;

      if ( bcf_update_genotypes(bcfHeader, bcfRecord, gts, p->M) ) die("Could not update GT field\n");

      // example of adding INFO fields
      bcf_update_info_int32(bcfHeader, bcfRecord, "AC", &ac[1], 1) ;
      bcf_update_info_int32(bcfHeader, bcfRecord, "AN", &an, 1) ;

      //write and progress
      bcf_write(fp, bcfHeader, bcfRecord) ;
      bcf_clear(bcfRecord) ;

      pbwtCursorForwardsRead(u) ;
    }

  // cleanup
  free(hap) ;
  free(gts) ;
  /*  pbwtCursorDestroy(u) ; */
  bcf_hdr_destroy(bcfHeader) ;
  bcf_destroy1(bcfRecord);
  hts_close(fp) ;
}

/******* end of file ********/

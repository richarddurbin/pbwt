/*  File: pbwtHtslib.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013
 *-------------------------------------------------------------------
 * Description: all the pbwt stuff that uses htslib, e.g. reading/writing vcf or bcf files
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  7 16:25 2015 (rd)
 * * Sep 22 23:03 2014 (rd): change for 64bit arrays
 * Created: Thu Oct 17 12:20:04 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>
#include <ctype.h>		/* for toupper() */

const char *pbwtHtslibVersionString(void)
{
    return hts_version();
}

static void readVcfSamples (PBWT *p, bcf_hdr_t *hr)
{
  int i, k ;

  p->samples = arrayCreate (p->M, int) ;
  for (i = 0 ; i < p->M/2 ; ++i)
    { k = sampleAdd (hr->samples[i]) ;
      array(p->samples, 2*i, int) = k ; /* assume diploid - could be cleverer */
      array(p->samples, 2*i+1, int) = k ;
    }
}

static int variation (const char *ref, const char *alt)
{
  int var, len = strlen (ref) + strlen (alt) + 2;
  char *buf = (char*) malloc(len);
  sprintf (buf, "%s\t%s", ref, alt) ;
  char *ptr = buf;
  while ( *ptr ) { *ptr = toupper(*ptr); ptr++; }
  dictAdd (variationDict, buf, &var) ;
  free(buf);
  return var ;
}


PBWT *pbwtReadVcfGT (char *filename, int isXY)  /* read GTs from vcf/bcf using htslib */
{
  int i, j, nHaplotypes, nSamplesKeep ;
  PBWT *p ;

  bcf_srs_t *sr = bcf_sr_init() ;
  if (!bcf_sr_add_reader (sr, filename)) die ("failed to open good vcf file\n") ;

  bcf_hdr_t *hr = sr->readers[0].header ;
  int nSamples = bcf_hdr_nsamples(hr) ;

  Array samples = arrayCreate (nSamples, int) ;
  Array ploidy = arrayCreate (nSamples, int) ;
  
  // read in samples from the header
  // use the isX/isY flags plus sample
  // sex information loaded previously
  // to determine number of haplotypes
  // and assign a per-sample ploidy
  for (i = 0, nHaplotypes = 0, nSamplesKeep = 0 ; i < nSamples ; ++i)
    {
      int k = array(samples,i,int) = sampleAdd (hr->samples[i]) ;
      Sample *s = sample (k) ; 
      array(ploidy, i, int) = 0 ;
      if (isXY==1 && s->isFemale) continue ;
      array(ploidy, i, int)++ ;
      nSamplesKeep++ ;
      nHaplotypes++ ;
      if (isXY==1 || (isXY==2 && s->isMale)) continue ;
      nHaplotypes++ ;
      array(ploidy, i, int)++ ;
    }
  // create the PBWT
  p = pbwtCreate(nHaplotypes, 0) ;
  p->samples = arrayReCreate(p->samples, p->M, int) ;
  if (isXY==2) p->isX = TRUE ;
  if (isXY==1) p->isY = TRUE ;

  // fill in the p->samples array
  for (i = 0, j = 0 ; i < nSamples ; ++i)
    {
      int iploidy = arr(ploidy, i, int) ;
      if (!iploidy) continue ; // sample ploidy 0, do not store in PBWT
      array(p->samples, j++, int) = arr(samples,i,int) ;
      if (iploidy==1) continue ;
      array(p->samples, j++, int) = arr(samples,i,int) ;
    }
  arrayDestroy (samples) ;

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

      // get a copy of GTs
      int ngt = bcf_get_genotypes(hr, line, &gt_arr, &mgt_arr) ;
      if (ngt <= 0) continue ;  // it seems that -1 is used if GT is not in the FORMAT
      if (ngt != nSamples && ngt != 2*nSamples) die ("%d != %d GT values at %s:%d - not haploid or diploid?",  ngt, nSamples, chrom, pos) ;
      
      memset (xMissing, 0, p->M) ;
      long wasMissing = nMissing ;

      /* copy the genotypes into array x[] */
      if (ngt == nSamples) 
        {
          /*
            All GTs are stored as haploid:
            - sample ploidy 0, do not store in PBWT
            - sample ploidy 2, treat haploid genotype as diploid homozygous A/A
            - sample ploidy 1, treat haploid genotype as is
          */
          for (i = 0, j = 0 ; i < nSamples ; i++)
            {
              int iploidy = arr(ploidy, i, int) ;
              if (!iploidy) continue ; // sample ploidy 0, do not store in PBWT
              if (iploidy==2) // sample ploidy 2, treat haploid genotype as diploid homozygous A/A
                {
                  if (gt_arr[i] == bcf_gt_missing)
                    { x[j] = 0 ;
                      x[j+1] = 0; /* use ref for now */
                      xMissing[j] = 1 ;
                      xMissing[j+1] = 1;
                      nMissing+=2 ;
                    }
                  else {
                    x[j] = bcf_gt_allele(gt_arr[i]) ;  // convert from BCF binary to 0 or 1
                    x[j+1] = x[j] ;  // convert from BCF binary to 0 or 1
                  }
                }
              else // sample ploidy 1, treat haploid genotype as is
                {
                  if (gt_arr[i] == bcf_gt_missing)
                    { x[j] = 0 ; /* use ref for now */
                      xMissing[j] = 1 ;
                      ++nMissing ;
                    }
                  else 
                    x[j] = bcf_gt_allele(gt_arr[i]) ;  // convert from BCF binary to 0 or 1
                }
              j += iploidy ;
            }
        }
      else
        {
          /*
            GTs are stored as diploid or mixture of diploid/haploid:
            - sample ploidy 0, do not store in PBWT
            - sample ploidy 2 and GT diploid, treat as is
            - sample ploidy 2 and GT haploid, convert to diploid homozygous A/A
            - sample ploidy 1 and GT haploid, treat as is
            - sample ploidy 1 and GT diploid homozygous, set to haploid A
            - sample ploidy 1 and GT diploid heterozygous, set to missing
          */
          for (i = 0, j = 0 ; i < ngt ; i+=2)
            {
              int iploidy = arr(ploidy, i/2, int) ;
              if (!iploidy) continue ; // sample ploidy 0, do not store in PBWT
              if (iploidy==2)
                { if (gt_arr[i] != bcf_gt_missing && gt_arr[i+1] == bcf_int32_vector_end) 
                    {
                      x[j]   = bcf_gt_allele(gt_arr[i]); // sample ploidy 2 and GT haploid, convert to diploid homozygous A/A
                      x[j+1] = bcf_gt_allele(gt_arr[i]);
                    }
                  else if (gt_arr[i] == bcf_gt_missing)
                    {
                      x[j] = 0 ; /* use ref for now */
                      x[j+1] = 0 ;
                      xMissing[j] = 1 ;
                      xMissing[j+1] = 1 ;
                      ++nMissing ;
                    }
                  else
                  {
                    x[j] = bcf_gt_allele(gt_arr[i]) ;  // convert from BCF binary to 0 or 1
                    x[j+1] = bcf_gt_allele(gt_arr[i+1]) ;  // convert from BCF binary to 0 or 1
                  }
                }
              else // if sample ploidy marked as haploid, treat haploid genotypes as is
                {
                  int a1 = bcf_gt_allele(gt_arr[i]) ;
                  int a2 = bcf_gt_allele(gt_arr[i+1]) ;
                  if (a1 == a2 || gt_arr[i+1] == bcf_int32_vector_end) 
                    {
                      if (gt_arr[i] == bcf_gt_missing)
                        {
                          x[j] = 0 ; /* use ref for now */
                          xMissing[j] = 1 ;
                          ++nMissing ;
                        }
                      else
                        x[j] = a1 ;  // convert from BCF binary to 0 or 1
                    }
                  else
                  {
                    x[j] = a1 ;
                    xMissing[j] = 1 ;
                    ++nMissing ;
                  }
                }
              j += iploidy ;
            }
        }

      BOOL no_alt = line->n_allele == 1;
      int n_allele = no_alt ? 2 : line->n_allele;

      /* split into biallelic sites filling in as REF ALT alleles */
      /* not in the REF/ALT site */
      for (i = 1 ; i < n_allele ; i++)
        {
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
          s->varD = variation (line->d.allele[0], no_alt ? "." : line->d.allele[i]) ;
        }

      if (nCheckPoint && !(p->N % nCheckPoint))  pbwtCheckPoint (u, p) ;
    }


  fflush(stdout);

  pbwtCursorToAFend (u, p) ;

  arrayDestroy (ploidy) ;
  if (gt_arr) free (gt_arr) ;
  bcf_sr_destroy (sr) ;
  free (x) ;
  pbwtCursorDestroy (u) ;  
  free (xMissing) ;

  fprintf (logFile, "read genotypes from %s with %d sample names and %ld sites on chromosome %s: M, N are %d, %d\n", 
         filename, nSamplesKeep, arrayMax(p->sites), p->chrom, p->M, p->N) ;
  if (p->missingOffset) fprintf (logFile, "%ld missing values at %d sites\n", 
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

  if (!p) die ("pbwtWriteVcf called without a valid pbwt") ;
  if (!p->sites) die ("pbwtWriteVcf called without sites") ;
  if (!p->samples) fprintf (logFile, "Warning: pbwtWriteVcf called without samples... using fake sample names PBWT0, PBWT1 etc...\n") ;
  fp = hts_open(filename,mode) ;
  if (!fp) die ("could not open file for writing: %s", filename) ;
  BOOL isDosage = p->dosageOffset ? TRUE : FALSE ;

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
  ksprintf(&str, "##pbwtVersion=%d.%d%s%s+htslib-%s", pbwtMajorVersion, pbwtMinorVersion, 
          strcmp(pbwtCommitHash(),"")==0 ? "" : "-", pbwtCommitHash(), pbwtHtslibVersionString()) ;
  bcf_hdr_append(bcfHeader, str.s) ;
  str.l = 0;
  ksprintf(&str, "##pbwtCommand=pbwt %s", commandLine) ;
  bcf_hdr_append(bcfHeader, str.s) ;
  free(str.s) ;
  bcf_hdr_append(bcfHeader, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">") ;
  bcf_hdr_append(bcfHeader, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">") ;
  bcf_hdr_append(bcfHeader, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">") ;
  if (isDosage)
    {
      bcf_hdr_append(bcfHeader, "##INFO=<ID=RefPanelAF,Number=A,Type=Float,Description=\"Allele frequency in imputation reference panel\">") ;
      bcf_hdr_append(bcfHeader, "##INFO=<ID=DR2,Number=A,Type=Float,Description=\"Estimated haploid dosage r^2 from imputation\">") ;
      bcf_hdr_append(bcfHeader, "##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Site was genotyped prior to imputation\">") ;
      bcf_hdr_append(bcfHeader, "##FORMAT=<ID=ADS,Number=.,Type=Float,Description=\"Allele dosage per haplotype\">") ;
      bcf_hdr_append(bcfHeader, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">") ;
      bcf_hdr_append(bcfHeader, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">") ;
    }
  
  int i, j, k ;
  for (i = 0, j = 0 ; i < p->M ; i += pbwtSamplePloidy (p, i))
    {
      if (p->samples)
        bcf_hdr_add_sample(bcfHeader, sampleName(pbwtSample (p, i))) ;
      else
        {
          kstring_t sname = {0,0,0} ;
          ksprintf(&sname, "PBWT%d", j++) ;
          bcf_hdr_add_sample(bcfHeader, sname.s) ;
          free(sname.s) ;
        }
    }
  bcf_hdr_add_sample(bcfHeader, 0) ; /* required to update internal structures */
  bcf_hdr_write(fp, bcfHeader) ;
  int nSamples = bcf_hdr_nsamples(bcfHeader) ;
  BOOL allHaploid = nSamples == p->M ; // all haploid, nSamples == nHaplotypes
  int nG = allHaploid ? 2*nSamples : 3*nSamples ; // length of a Number=G FORMAT vector

  bcf1_t *bcfRecord = bcf_init1() ;
  uchar *hap = myalloc (p->M, uchar) ;
  int32_t *gts = myalloc (allHaploid ? nSamples : 2*nSamples, int32_t) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  double *d = 0 ;
  double *gps = isDosage ? myalloc (nG, double) : NULL;
  double *ds = isDosage ? myalloc (nSamples, double) : NULL;
  double *ad = isDosage ? myalloc (p->M, double) : NULL;
  float *fls = myalloc (3*nSamples, float);
  uchar *missing = mycalloc (p->M, uchar) ;
  if (!p->missingOffset) bzero (missing, p->M) ;

  for (i = 0 ; i < p->N ; ++i)
    {
      Site *s = arrp(p->sites, i, Site) ;
      bcf_float_set_missing(bcfRecord->qual) ;
      bcfRecord->rid = bcf_hdr_name2id(bcfHeader, p->chrom) ;
      bcfRecord->pos = s->x - 1 ;
      char *als = strdup( dictName(variationDict, s->varD) ), *ss = als ;
      while ( *ss ) { if ( *ss=='\t' ) *ss = ',' ; ss++ ; }
      bcf_update_alleles_str(bcfHeader, bcfRecord, als) ;
      free(als) ;
      bcf_add_filter(bcfHeader, bcfRecord, bcf_hdr_id2int(bcfHeader, BCF_DT_ID, "PASS")) ;

      if (p->missingOffset)
        {
          if (!arr(p->missingOffset, i, long)) bzero (missing, p->M) ;
          else unpack3 (arrp(p->zMissing, arr(p->missingOffset,i,long), uchar), p->M, missing, 0) ;
        }

      if (isDosage) d = pbwtDosageRetrieve (p, u, d, i) ;
      // map haplotypes and dosages to sample order
      for (j = 0 ; j < p->M ; ++j)
        {
          hap[u->a[j]] = u->y[j] ;
          if (isDosage) ad[u->a[j]] = d[j] ;
        }

      /*
        The dosages d are actually posterior probabilities per haplotype
        Diploid genotype dosage:
          gd[n] = d[2*n] + d[2*n+1]
        Haploid genotype dosage:
          gd[n] = d[n]
        Diploid genotype posterior probabilities:
          gp[n][0] = (1-d[2*n]) * (1-d[2*n+1])
          gp[n][1] = d[2*n] + d[2*n+1] - 2*d[2*n]*d[2*n+1]
          gp[n][2] = d[2*n] * d[2*n+1]
        Haploid genotype posterior probabilities:
          gp[n][0] = 1-d[n]
          gp[n][1] = d[n]
      */

      int ac[2] = {0,0};
      float raf = s->refFreq;
      float info = s->imputeInfo;
      if (allHaploid)
        {
          for (j = 0 ; j < p->M ; j++)
            {
              if (isDosage)
                {
                   ds[j] = ad[j] ;
                   gps[2*j] = 1-ad[j] ;
                   gps[2*j+1] = ad[j] ;
                }
              if (hap[j]<0)
                gts[j] = bcf_gt_missing ;
              else
                {
                  gts[j] = bcf_gt_unphased(hap[j]) ;
                  ac[hap[j]]++ ;
                }
            }
        }
      else
        {
          for (j = 0, k = 0 ; j < p->M ; k+=2)
            {
              int ploidy = pbwtSamplePloidy (p, j) ;
              if (isDosage)
                {
                  ds[k/2] = ploidy==1 ? ad[j] : ad[j] + ad[j+1] ;
                  gps[3*k/2] = ploidy==1 ? 1-ad[j] : (1-ad[j]) * (1-ad[j+1]) ;
                  gps[3*k/2+1] = ploidy==1 ? ad[j] : ad[j] + ad[j+1] - 2*ad[j]*ad[j+1] ;
                  gps[3*k/2+2] = ploidy==1 ? -1 : ad[j] * ad[j+1] ; // encode vector end as -1 for conversion to bcf_float_vector_end later
                }
              if (missing[j])
                gts[k] = bcf_gt_missing ;
              else
                {
                  gts[k] = bcf_gt_unphased(hap[j]) ;
                  ac[hap[j]]++ ;
                }
              if (ploidy==2)
                {
                  if (missing[j+1])
                    gts[k+1] = bcf_gt_missing ;
                  else
                    {
                      gts[k+1] = p->isUnphased ? bcf_gt_unphased(hap[j+1]) : bcf_gt_phased(hap[j+1]) ;
                      ac[hap[j+1]]++ ;
                    }
                }
              else
                  gts[k+1] = bcf_int32_vector_end ;

              j += ploidy ;
            }
        }
      int an = ac[0] + ac[1] ;

      if ( bcf_update_genotypes(bcfHeader, bcfRecord, gts, allHaploid ? nSamples : 2*nSamples) ) die("Could not update GT field\n");
      if (p->isRefFreq)
          if ( bcf_update_info_float(bcfHeader, bcfRecord, "RefPanelAF", &raf, 1) ) die("Could not update INFO/RefPanelAF field\n") ;
      if (isDosage)
        {
          if ( bcf_update_info_float(bcfHeader, bcfRecord, "DR2", &info, 1) ) die("Could not update INFO/DS field\n") ;
          int k ;
          // dosages stored as double, but BCF required floats
          // may be a better way to handle this, but it works
          // for the moment
          if (allHaploid)
            {
              for (k = 0 ; k < nSamples ; ++k)
                fls[k] = (float)(ad[k]) ;
              if ( bcf_update_format_float(bcfHeader, bcfRecord, "ADS", fls, nSamples) ) die("Could not update FORMAT/ADS field\n") ;
            }
          else
          {
            for (j = 0, k = 0 ; j < p->M ; k+=2)
              {
                int ploidy = pbwtSamplePloidy (p, j) ;
                fls[k] = (float)(ad[j++]) ;
                if (ploidy==1)
                  bcf_float_set_vector_end(fls[k+1]) ;
                else
                  fls[k+1] = (float)(ad[j++]) ;
              }
            if ( bcf_update_format_float(bcfHeader, bcfRecord, "ADS", fls, 2*nSamples) ) die("Could not update FORMAT/ADS field\n") ;
          }
          for (j = 0 ; j < nSamples ; ++j)
            fls[j] = (float)(ds[j]) ;
          if ( bcf_update_format_float(bcfHeader, bcfRecord, "DS", fls, nSamples) ) die("Could not update FORMAT/DS field\n") ;
          for (j = 0 ; j < nG ; ++j) {
            if (gps[j]<0)
              bcf_float_set_vector_end(fls[j]) ;
            else
              fls[j] = (float)(gps[j]) ;
          }
          if ( bcf_update_format_float(bcfHeader, bcfRecord, "GP", fls, nG) ) die("Could not update FORMAT/GP field\n") ;
        }

      bcf_update_info_int32(bcfHeader, bcfRecord, "AC", &ac[1], 1) ;
      bcf_update_info_int32(bcfHeader, bcfRecord, "AN", &an, 1) ;
      if (!s->isImputed) bcf_update_info_flag(bcfHeader, bcfRecord, "TYPED", NULL, 1) ;

      //write and progress
      bcf_write(fp, bcfHeader, bcfRecord) ;
      bcf_clear(bcfRecord) ;

      pbwtCursorForwardsRead(u) ;
    }

  // cleanup
  free(d) ;
  free(hap) ;
  free(gts) ;
  free(missing) ;
  if (isDosage) { free(fls) ; free(gps) ; free(ds) ; free(ad) ; }
  pbwtCursorDestroy(u) ; // this was commented, not sure why??
  bcf_hdr_destroy(bcfHeader) ;
  bcf_destroy1(bcfRecord);
  hts_close(fp) ;

  fprintf (logFile, "written vcf file: %d records, %d samples, and %d haplotypes\n", p->N, nSamples, p->M) ;
}

/******* end of file ********/

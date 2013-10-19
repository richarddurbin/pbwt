/*  File: pbwtHtslib.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013
 *-------------------------------------------------------------------
 * Description: all the pbwt stuff that uses htslib, e.g. reading/writing vcf or bcf files
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 19 14:52 2013 (rd)
 * Created: Thu Oct 17 12:20:04 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "pbwt.h"
#include <htslib/synced_bcf_reader.h>

PBWT *pbwtReadVcf (char *filename)	/* read vcf/bcf using htslib */
{
  PBWT *p ;
  int i, j ;

  bcf_srs_t *sr = bcf_sr_init() ;
  bcf_sr_add_reader (sr, filename) ;

  bcf_hdr_t *hr = sr->readers[0].header ;
  p = pbwtCreate (bcf_hdr_nsamples(hr)*2) ; /* assume diploid! */
  for (i = 0 ; i < p->M/2 ; ++i)
    /* add sample hr->samples[i] */ ;

  int mgt_arr = 0, *gt_arr = NULL;
  int mpl_arr = 0, *pl_arr = NULL;
  while (bcf_sr_next_line (sr)) 
    { bcf1_t *line = bcf_sr_get_line(sr,0) ;
      const char* chrom = bcf_seqname(hr,line) ;
      int pos = line->pos + 1 ;                      // coordinates are 0-based
      if ( line->n_allele !=2 ) continue;       // not a biallelic site
      const char *ref = line->d.allele[0] ;
      const char *alt = line->d.allele[1] ;

      printf("%s:%d %s %s", chrom,pos,ref,alt);

      // get a copy of GTs
      int ngt = bcf_get_genotypes(hr, line, &gt_arr, &mgt_arr) ;
      if (!ngt) continue ;             // GT not present
      ngt /= bcf_hdr_nsamples (hr) ;      // ngt is now the number of values per sample
      for (i = 0 ; i < bcf_hdr_nsamples (hr) ; i++)
      {
        // skip missing genotypes and haploid samples
        if ( gt_arr[i*ngt]==bcf_gt_missing ) continue;  // missing
        if ( gt_arr[i*ngt+1]==bcf_gt_missing || gt_arr[i*ngt+1]==bcf_int32_vector_end ) continue; // missing or haploid
        int al1 = bcf_gt_allele(gt_arr[i*ngt]);         // convert from BCF binary representation to 0 or 1
        int al2 = bcf_gt_allele(gt_arr[i*ngt+1]);       //      "  "
        assert( al1==0 || al1==1 );
        assert( al2==0 || al2==1 );

        printf(" %d/%d", al1,al2);
      }

      // get a copy of the PL vectors
      int npl = bcf_get_format_int(hr, line, "PL", &pl_arr, &mpl_arr);
      if ( !npl ) continue;             // PL not present
      npl /= bcf_hdr_nsamples(hr);     // number of values per samples
      for (i=0; i<bcf_hdr_nsamples(hr); i++)    // iterate over samples
      {
        for (j=0; j<npl; j++)         // iterate over PL values (genotypes)
        {
          // check for shorter vectors (haploid genotypes amongst diploid genotypes)
          if ( pl_arr[i*npl+j]==bcf_int32_vector_end ) break;

          // skip missing values
          if ( pl_arr[i*npl+j]==bcf_int32_missing ) continue;

          // do something with j-th PL value
          printf(" %d", pl_arr[i*npl+j]);
        }
      }

      printf("\n");
    }
  free(gt_arr);
  free(pl_arr);
  bcf_sr_destroy(sr);

  return p ;
}

/******* end of file ********/

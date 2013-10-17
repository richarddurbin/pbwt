/*  File: pbwtHtslib.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2013
 *-------------------------------------------------------------------
 * Description: all the pbwt stuff that uses htslib, e.g. reading/writing vcf or bcf files
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 17 12:51 2013 (rd)
 * Created: Thu Oct 17 12:20:04 2013 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "pbwt.h"
#include <htslib/synced_bcf_reader.h>

PBWT *pbwtReadVcf (char *filename)	/* read vcf/bcf using htslib */
{
  PBWT *p ;
  int i ;

  bcf_srs_t *sr = bcf_sr_init() ;
  bcf_sr_add_reader (sr, filename) ;

  bcf_hdr_t *hr = sr->readers[0].header ;
  p = pbwtCreate (bcf_hdr_nsamples(hr)*2) ; /* assume diploid! */
  for (i = 0 ; i < p->M/2 ; ++i)
    /* add sample hr->samples[i] */ ;

  while (bcf_sr_next_line (sr)) 
    { bcf1_t *line = bcf_sr_get_line(sr,i) ;
      char* chrom = bcf_seqname(hr,line) ;
    }

  return p ;
}

/******* end of file ********/

/*  File: pbwtPaint.c
 *  Author: Richard Durbin (rd@sanger.ac.uk) and  Daniel Lawson (dan.lawson@bristol.ac.uk)
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
 * Description: tools for chromosome painting as in ChromoPainter, FineStructure etc.
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  4 23:55 2017 (rd)
 * Created: Tue Apr  1 11:34:41 2014 (rd)
 *-------------------------------------------------------------------
 */

//#define LEGACY_MATRIX

#include "zlib.h"
#include "pbwt.h"
//#include "dynamiclist.h"

typedef struct { int j ; int start ; int end ; } MatchSegment ;
/* matches are semi-open [start,end) so length is end-start */
static Array *maxMatch = 0 ;	/* arrays of MatchSegment */
static void reportMatch (int i, int j, int start, int end)
{ MatchSegment *ms = arrayp (maxMatch[i], arrayMax(maxMatch[i]), MatchSegment) ;
  ms->j = j ; ms->start = start ; ms->end = end ;
}

static inline int isprevind(int i,int *map_indhap){
  if(i==0) return(0);
  return(map_indhap[i]==map_indhap[i-1]);
}

static inline void printAll(int ii,int Ninds,
			    double *t_counts,double *t_counts2,double *t_counts3,double *t_totlengths,double nregions,
			    gzFile fc,gzFile fc2,gzFile fc3,gzFile fl,gzFile fr){
  int jj ; for (jj = 0 ; jj < Ninds ; ++jj) {
    if(t_counts[jj]){
      gzprintf (fc, "IND%i IND%i %.4f\n", ii+1,jj+1,t_counts[jj]) ; 
      gzprintf (fl, "IND%i IND%i %.4f\n", ii+1,jj+1,t_totlengths[jj]) ; 
      gzprintf (fc2,"IND%i IND%i %.4f\n", ii+1,jj+1,t_counts2[jj]) ; 
      gzprintf (fc3,"IND%i IND%i %.4f\n", ii+1,jj+1,t_counts3[jj]) ; 
    }
  }
  gzprintf (fr,"IND%i %.2f\n",ii+1, nregions) ; 
}

void paintAncestryMatrix (PBWT *p, char* fileRoot,int chunksperregion)
{
  int ploidy=2;
  int Ninds=p->M/ploidy;
  double **totlengths = 0 ;
  double **counts = 0 ;
  double **counts2 = 0 ; /* store sums of squares of counts in 100 block bins */
  double **counts3 = 0 ; /* store sums of counts in 100 block bins */
  double *nregions = 0; /* store number of regions found */

  int i, j, k ;
  int *map_indhap = mycalloc (p->M,int);
  for (i = 0 ; i < p->M ; ++i) map_indhap[i]=i/ploidy;

  double *totCounts = mycalloc (Ninds, double) ;

  nregions = mycalloc (Ninds,double);
  counts = myalloc (Ninds, double*) ; counts2 = myalloc (Ninds, double*) ; counts3 = myalloc (Ninds, double*) ; totlengths= myalloc (Ninds, double*);
  for (i = 0 ; i < Ninds ; ++i) 
    { 
      counts[i] = mycalloc (Ninds, double) ;
      counts2[i] = mycalloc (Ninds, double) ;
      counts3[i] = mycalloc (Ninds, double) ;
      totlengths[i] = mycalloc (Ninds, double) ;
    }
  memset (nregions, 0, sizeof(double)*Ninds) ;

  maxMatch = myalloc (p->M, Array) ;
  for (i = 0 ; i < p->M ; ++i) maxMatch[i] = arrayCreate (1024, MatchSegment) ;
  matchMaximalWithin (p, reportMatch) ;  /* store maximal matches in maxMatch */
  double *partCounts = myalloc (Ninds, double) ;
  /* now weight per site based on distance from ends */

  for (i = 0 ; i < p->M ; ++i)
    { 
      //      printf("Processing individual %i (haplotype %i)\n",i/ploidy,i);
      MatchSegment *m1 = arrp(maxMatch[i],0,MatchSegment), *m ;
      int n1 = 1 ;		/* so don't have an empty chunk to start with! */
      MatchSegment *mStop = arrp(maxMatch[i], arrayMax(maxMatch[i])-1, MatchSegment) ;
      memset (partCounts, 0, sizeof(double)*Ninds) ;
      for (k = 1 ; k < p->N ; k++)
	{ double sum = 0 ;
	  while (m1->end <= k && m1 < mStop)
	    { if ((n1 % chunksperregion)==0)
		{ int jj ; for (jj = 0 ; jj < Ninds ; ++jj) {
		    counts2[map_indhap[i]][jj] += partCounts[jj]*partCounts[jj] ;
		    counts3[map_indhap[i]][jj] += partCounts[jj] ;
		  }
		  memset (partCounts, 0, sizeof(double)*Ninds) ;
		  nregions[map_indhap[i]]+=1.0;
		}
	      ++m1 ; ++n1 ;
	    }
	  for (m = m1 ; m->start < k && m <= mStop ; ++m) 
	    sum += (k - m->start) * (m->end - k) ;
	  if (sum)
	    for (m = m1 ; m->start < k && m <= mStop ; ++m) {
	      totlengths[map_indhap[i]][map_indhap[m->j]] += (k - m->start) * (m->end - k) / sum;
 	      double thiscount=(k - m->start) * (m->end - k) / sum/(m->end - m->start);
	      counts[map_indhap[i]][map_indhap[m->j]] += thiscount;
	      partCounts[map_indhap[m->j]] += thiscount;
	    }
	}
    }
  
  free(map_indhap);
  free (partCounts) ;

  /* report results */
  FILE *fc = fopenTag (fileRoot, "chunkcounts.out", "w") ;
  FILE *fl = fopenTag (fileRoot, "chunklengths.out", "w") ;
  FILE *fc2 = fopenTag (fileRoot, "regionsquaredchunkcounts.out", "w") ;
  FILE *fc3 = fopenTag (fileRoot, "regionchunkcounts.out", "w") ;
  fprintf (fc,"RECIPIENT") ; 
  fprintf (fl,"RECIPIENT") ; 
  fprintf (fc2,"RECIPIENT nregions") ; 
  fprintf (fc3,"RECIPIENT nregions") ; 
  for (i = 0 ; i < Ninds ; ++i)    {
    fprintf (fc," IND%i",i+1) ; 
    fprintf (fl," IND%i",i+1) ; 
    fprintf (fc2," IND%i",i+1) ; 
    fprintf (fc3," IND%i",i+1) ; 
  }
  fputc ('\n', fc) ;
  fputc ('\n', fl) ;
  fputc ('\n', fc2) ;
  fputc ('\n', fc3) ;
 
 for (i = 0 ; i < Ninds ; ++i)    {
   fprintf (fc3,"IND%i %.2f",i+1, nregions[i]) ; 
   fprintf (fc2,"IND%i %.2f",i+1, nregions[i]) ; 
   fprintf (fl,"IND%i",i+1) ; 
   fprintf (fc,"IND%i",i+1) ; 
 for (j = 0 ; j < Ninds ; ++j) 
	{ 
	  fprintf (fc, " %.4f", counts[i][j]) ; 
	  fprintf (fl, " %.4f", totlengths[i][j]) ; 
 	  fprintf (fc2," %.4f", counts2[i][j]) ; 
	  fprintf (fc3," %.4f", counts3[i][j]) ; 
	  totCounts[i] += counts[i][j] ; 
	}
      fputc ('\n', fc) ;
      fputc ('\n', fl) ;
      fputc ('\n', fc2) ;
      fputc ('\n', fc3) ;
      if (isCheck && (i%2) && p->samples) 
	fprintf (logFile, "%s %8.4g %8.4g\n", 
		 sampleName (sample(p,i-1)), totCounts[i-1], totCounts[i]) ;
    }
  fclose (fc) ; fclose (fl) ; fclose (fc2) ;fclose (fc3) ;
  timeUpdate(logFile);
  /* clean up */
  for (i = 0 ; i < Ninds ; ++i) { free (counts[i]) ; free (counts2[i]) ; free (counts3[i]) ; free (totlengths[i]) ; }
  free (counts) ; free (counts2) ; free (counts3) ; free (totCounts) ; free (nregions); free(totlengths);
}

void paintAncestryMatrixSparse (PBWT *p, char* fileRoot,int chunksperregion,int cutoff)
{
  int i, j, k ;
  int ploidy=2;
  int Ninds=p->M/ploidy;
  int *map_indhap = mycalloc (p->M,int);
  for (i = 0 ; i < p->M ; ++i) map_indhap[i]=i/ploidy;

  double *nregions = 0; /* store number of regions found */
  nregions = mycalloc (Ninds,double);
  //  double *totCounts = mycalloc (Ninds, double) ;

  gzFile fr = gzopenTag (fileRoot, "nregions.s.out.gz", "w") ;
  gzFile fc = gzopenTag (fileRoot, "chunkcounts.s.out.gz", "w") ;
  gzFile fl = gzopenTag (fileRoot, "chunklengths.s.out.gz", "w") ;
  gzFile fc2 = gzopenTag (fileRoot, "regionsquaredchunkcounts.s.out.gz", "w") ;
  gzFile fc3 = gzopenTag (fileRoot, "regionchunkcounts.s.out.gz", "w") ;
  /// The dynamic structure storing the sparse matrices

  int *t_obs = myalloc (Ninds,int);
  double *t_counts = myalloc (Ninds,double);
  double *t_counts2 = myalloc (Ninds,double);
  double *t_counts3 = myalloc (Ninds,double);
  double *t_totlengths = myalloc (Ninds,double);

  maxMatch = myalloc (p->M, Array) ;
  for (i = 0 ; i < p->M ; ++i) maxMatch[i] = arrayCreate (1024, MatchSegment) ;
  matchMaximalWithin (p, reportMatch) ;  /* store maximal matches in maxMatch */

  double *partCounts = myalloc (Ninds, double) ; 

  /* now weight per site based on distance from ends */
  /// i =1 .. M is each SNP
  for (i = 0 ; i < p->M ; ++i)
    {
      //      printf("Processing Individual %i in haplotype %i\n",i/ploidy,i);
      MatchSegment *m1 = arrp(maxMatch[i],0,MatchSegment);      // m1 is a SEGMENT that matches at SNP i
      MatchSegment *m ;// m is another segment
      int n1 = 1 ; // number of chunks found so far. 1 so don't have an empty chunk to start with! 

      MatchSegment *mStop = arrp(maxMatch[i], arrayMax(maxMatch[i])-1, MatchSegment) ;//
      
      // Clear records if we have a new individual
      if(!isprevind(i,map_indhap)){
		if(i>0) printAll(map_indhap[i-1],Ninds,t_counts,t_counts2,t_counts3,t_totlengths,nregions[map_indhap[i-1]],fc,fc2,fc3,fl,fr);
	memset (t_obs, 0, sizeof(int)*Ninds) ;
	memset (partCounts, 0, sizeof(double)*Ninds) ;
	memset (t_counts, 0, sizeof(double)*Ninds) ;
	memset (t_counts2, 0, sizeof(double)*Ninds) ;
	memset (t_counts3, 0, sizeof(double)*Ninds) ;
	memset (t_totlengths, 0, sizeof(double)*Ninds) ;
      }

      for (k = 1 ; k < p->N ; k++) 
	{ double sum = 0 ;
	  while (m1->end <= k && m1 < mStop)
	    { if ((n1 % chunksperregion)==0) { 
		int jj ; for (jj = 0 ; jj < Ninds ; ++jj) {
		  if(partCounts[jj]){
		    t_counts2[jj] += partCounts[jj]*partCounts[jj] ;
		    t_counts3[jj] += partCounts[jj] ;
		  }
		}
		memset (partCounts, 0, sizeof(double)*Ninds) ;
		nregions[map_indhap[i]]+=1.0;
	      }
	      ++m1 ; ++n1 ;
	    }
	  // for every individual who has a sufficiently long match
	  for (m = m1 ; m->start < k && m <= mStop ; ++m) 
	    sum += (k - m->start) * (m->end - k) ;
	  if (sum)
	    for (m = m1 ; m->start < k && m <= mStop ; ++m) {
	      double thislengths = (k - m->start) * (m->end - k) / sum;
 	      double thiscount=(k - m->start) * (m->end - k) / sum/(m->end - m->start);

	      t_obs[map_indhap[m->j]] = 1;
	      t_totlengths[map_indhap[m->j]] += thislengths;
	      t_counts[map_indhap[m->j]] += thiscount;
	      partCounts[map_indhap[m->j]] += thiscount;
	      //	      dynamicListSet(usedinds,map_indhap[m->j],1);
	    }
	} // end loop over the N loci
      
    }// end loop over individuals

    printAll(map_indhap[p->M-1],Ninds,
	   t_counts,t_counts2,t_counts3,t_totlengths,nregions[map_indhap[p->M-1]],
	   fc,fc2,fc3,fl,fr);  

  /* clean up */
  free (t_obs) ;
  free (t_counts) ;
  free (t_counts2) ;
  free (t_counts3) ;
  free (t_totlengths) ;

  free (partCounts) ;
  free (nregions) ;
  free(map_indhap);

   gzclose (fc) ; gzclose (fl) ;  gzclose (fc2) ;gzclose (fc3) ; gzclose(fr);
 


}


/* end of file */

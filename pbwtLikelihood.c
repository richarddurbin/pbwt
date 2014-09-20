/*  File: pbwtLikelihood.c
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
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 21 08:57 2014 (rd)
 * Created: Sat Apr 26 22:29:13 2014 (rd)
 *-------------------------------------------------------------------
 */

#include "pbwt.h"
#include <math.h>

/** line search code to maximise a function **/

double lineSearchPositive (double xInit, double tol, double (*function)(double))
/* find value to maximise function within tolerance */
{
  if (tol <= 1.0) die ("tolerance %f in lineSearchPostive() must be > 1.0", tol) ;
  double x0 = 0.9*xInit, y0 = (*function)(x0) ;
  double x1 = 1.1*xInit, y1 = (*function)(x1) ;
  double x2, y2 ;
  while (y0 < y1)
    { x2 = 3*x1 - 2*x0 ; if (x2 > 2.0*x1) x2 = 2.0*x1 ; y2 = (*function)(x2) ;
      if (isCheck) printf ("x0 %.4f %.4f > x1 %.4f %.4f\n", x0, y0, x1, y1) ;
      if (y1 > y2) break ;
      x0 = x1 ; y0 = y1 ; x1 = x2 ; y1 = y2 ;
    }
  while (y0 > y1)
    { if (isCheck) printf ("x0 %.4f %.4f < x1 %.4f %.4f\n", x0, y0, x1, y1) ;
      x2 = x1 ; y2 = y1 ; x1 = x0 ; y1 = y0 ;
      x0 = 3*x1 - 2*x2 ; if (x0 < 0.5*x1) x0 = 0.5*x1 ;
      y0 = (*function)(x0) ; 
    }
  /* now should have y1 > y0 and y1 > y2  */
  /* repeatedly fit a quadratic and pick the minimum */
  /* y = ax^2 - 2bx + c:  min at b/a */
  /* (y2-y1) = a(x2^2-x1^2) - 2b(x2-x1) */
  while (x2/x0 > tol)
    { double x ; 	/* new test value */
      if ((x1 - x0) > 2*(x2 - x1))
	{ x = 0.5*(x0 + x1) ; if (isCheck) printf ("split 01: ") ; }
      else if ((x2 - x1) > 2*(x1 - x0))
	{ x = 0.5*(x1 + x2) ; if (isCheck) printf ("split 12: ") ; }
      else
	{ double a = ((y2-y1)*(x1-x0) - (y1-y0)*(x2-x1))
	    / ((x2*x2-x1*x1)*(x1-x0) - (x1*x1-x0*x0)*(x2-x1)) ;
	  double b = 0.5 * (a * (x2*x2-x1*x1) - (y2-y1)) / (x2-x1) ;
	  x = b/a ;
	  if (isCheck) printf ("estimate: ") ;
	}
      double y = (*function)(x) ;
      if (isCheck) printf ("x/y0 %.4f %.4f  x/y1 %.4f %.4f x/y2 %.4f %.4f  x/ynew %.4f %.4f\n",
			   x0, y0, x1, y1, x2, y2, x, y) ;
      if (x > x1)
	if (y > y1) { x0 = x1 ; y0 = y1 ; x1 = x ; y1 = y ; }
	else { x2 = x ; y2 = y ; }
      else
	if (y > y1) { x2 = x1 ; y2 = y1 ; x1 = x ; y1 = y ; }
	else { x0 = x ; y0 = y ; }
    }
  return x1 ;
}

/*******************************************************************************/

static void simpleEntropy (PBWT *p)
{
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  int i, j, d ;
  double LL = 0, f ;
  long dTotStick = 0, nTotStick = 0 ;
  long dTotSwitch = 0, nTotSwitch = 0 ;

  for (i = 0 ; i < p->N ; ++i)
    { int last = u->y[0] ;
      for (j = 1 ; j < p->M ; ++j)
	{ d = i+1 - u->d[j] ;
	  if (u->y[j] == u->y[j-1]) { dTotStick += d ; ++nTotStick ; }
	  else { dTotSwitch += d ; ++nTotSwitch ; }
	}
      f = u->c/(double)p->M ;
      if (f > 0 && f < 1) LL += f * log(f) + (1-f) * log(1-f) ;
      pbwtCursorForwardsReadAD (u, i) ;
    }
  pbwtCursorDestroy (u) ;

  printf ("Fraction switch %.4f  av dStick %.1f av dSwitch %.1f\n", 
	  nTotSwitch / (double)(nTotStick+nTotSwitch),
	  dTotStick / (double)nTotStick, dTotSwitch / (double)nTotSwitch) ;
  
  printf ("Simple entropy per cell %f\n", LL/p->N) ;
}

/** package global variables we need for maximising likelihood for pbwt model **/

static Array info ;
static double alphaSearch, betaSearch ;

typedef struct {
  int nStick, nSwitch ;
} RowInfo ;

static Array buildRowInfo (PBWT *p, int MAX) /* array of RowInfo */
/* record how many times for each d we stick or switch going down the column */
{
  Array info = arrayCreate (4096, RowInfo) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  int i, j ;

  for (i = 0 ; i < p->N ; ++i)
    { for (j = 1 ; j < p->M ; ++j)
	{ int d = i+1 - u->d[j] ; if (d > MAX) d = MAX ;
	  if (u->y[j] == u->y[j-1]) 
	    ++arrayp(info,d,RowInfo)->nStick ; 
	  else 
	    ++arrayp(info,d,RowInfo)->nSwitch ;
	}
      pbwtCursorForwardsReadAD (u, i) ;
    }

  if (isStats)
    { int totStick = 0, totSwitch = 0 ; int lastStick = 0, lastSwitch = 0 ;
      double lastF = 1.0 ;
      for (i = 0 ; i < arrayMax(info) ; ++i)
	{ totStick += arrp(info,i,RowInfo)->nStick ;
	  totSwitch += arrp(info,i,RowInfo)->nSwitch ;
	  if (!((i+1)%100))
	    { double f = (totSwitch-lastSwitch)*100.0/(totStick+totSwitch-lastStick-lastSwitch) ;
	      printf ("%d  %d  %d  %.2f  %.3f\n", i+1, (totStick-lastStick), (totSwitch-lastSwitch), f, f/lastF) ; 
	      lastStick = totStick ; lastSwitch = totSwitch ; lastF = f ;
	    }
	}
      printf ("RowInfo counts: stick %d  switch %d", totStick, totSwitch) ; 
      printf (" %%stick %.2f  %%switch %.2f\n", 
	      totStick*100.0/(totStick+totSwitch), totSwitch*100.0/(totStick+totSwitch)) ;
    }

  pbwtCursorDestroy (u) ;

  return info ;
}

static double pbwtLogLikelihood (Array info, double alpha, double beta)
{
  int d ;
  double like = 0.0 ;
  RowInfo *inf = arrp(info,0,RowInfo) ;
  for (d = 0 ; d < arrayMax(info) ; ++d, ++inf) 
    { like += inf->nStick * log (1.0 - exp(-alpha*d - beta)) ;
      like += inf->nSwitch * (-alpha*d - beta) ;
    }

  return like ;
}

static double betaSearchLL (double beta)
{ return pbwtLogLikelihood (info, alphaSearch, beta) ; }

static double alphaSearchLL (double alpha)
{ alphaSearch = alpha ;
  betaSearch = lineSearchPositive (betaSearch, 1.001, betaSearchLL) ;
  return pbwtLogLikelihood (info, alphaSearch, betaSearch) ; 
}

/** drop one model: sum over sequences x_i LL(X)/LL(X\x_i) **/

typedef struct {
  int n[8] ;
  int nTot ;
} RowInfoDropOne ;

static Array buildRowInfoDropOne (PBWT *p, int MAX) /* array of RowInfoDropOne */
/* record how many times we see each set of 3 consecutive values in y, encoded as k,
   as a function of the pair of d values between them, encoded as dd */
{
  Array info = arrayCreate (4096, RowInfoDropOne) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  int i, j, k, d1, d2, dd ;

  for (i = 0 ; i < p->N ; ++i)
    { for (j = 0 ; j < p->M ; ++j)
	{ if (!u->d[j] || !u->d[j+1]) continue ; /* ignore edge effects */
	  if (j == 0) 
	    { k = (u->y[j] << 1) + u->y[j+1] ; d1 = 0 ; d2 = i+1 - u->d[j+1] ; }
	  else if (j < p->M-1)
	    { k = (u->y[j-1] << 2) + (u->y[j] << 1) + u->y[j+1] ; 
	      d1 = i+1 - u->d[j] ; d2 = i+1 - u->d[j+1] ; 
	    }
	  else
	    { k = (u->y[j-1] << 2) + (u->y[j] << 1) ; d1 = i+1 - u->d[j] ; d2 = 0 ; }
	  d1 /= 10 ; d2 /= 10 ;
	  if (d1 > MAX) d1 = MAX ; if (d2 > MAX) d2 = MAX ;
	  if (d1 < d2) dd = d2*d2 + d1 ; else dd = d1*d1 + d1 + d2 ;
	  arrayp(info,dd,RowInfoDropOne)->n[k] += 1 ; 
	  arrp(info,dd,RowInfoDropOne)->nTot += 1 ;
	}
      pbwtCursorForwardsReadAD (u, i) ;
    }

  if (isStats)
    { int kTot[8] ; for (k = 0 ; k < 8 ; ++k) kTot[k] = 0 ;
      for (dd = 0 ; dd < arrayMax(info) ; ++dd)
	if (arrp(info,dd,RowInfoDropOne)->nTot)
	  for (k = 0 ; k < 8 ; ++k) kTot[k] += arrp(info,dd,RowInfoDropOne)->n[k] ;
      printf ("RowInfoDropOne counts: ") ; 
      double tot = 0.0 ; 
      for (k = 0 ; k < 8 ; ++k) { printf (" %d", kTot[k]) ; tot += kTot[k] ; }
      printf (" %%stick %.1f  %%drift %.1f  %%flip %.1f\n", 
	      (kTot[0]+kTot[7])*100.0/tot, (kTot[1]+kTot[3]+kTot[4]+kTot[6])*100.0/tot, 
	      (kTot[2]+kTot[5])*100.0/tot) ;
    }

  pbwtCursorDestroy (u) ;

  return info ;
}

static double pbwtLogLikelihoodDropOne (Array info, double alpha, double beta)
/* this is the sum of "leave-one-out" likelihoods leaving each sequence out */
{
  int dmax = sqrt ((double)arrayMax(info)) ;
  double *pSwitch = myalloc (dmax+1, double) ; /* actually log_p values */
  double *pStick = myalloc (dmax+1, double) ; /* actually log_p values */
  int dd, d1, d2, dmin ; 
  for (d1 = 0 ; d1 <= dmax ; ++d1) 
    { pSwitch[d1] = - alpha*d1 - beta ; pStick[d1] = log (1.0 - exp(pSwitch[d1])) ; }
  double like = 0.0 ;
  RowInfoDropOne *inf = arrp(info,0,RowInfoDropOne) ;
  d1 = 0 ; d2 = 0 ;
  for (dd = 0 ; dd < arrayMax(info) ; ++dd, ++inf) 
    { if (inf->nTot)
	{ int *nn = inf->n ;
	  dmin = (d1 < d2) ? d1 : d2 ;
	  dmax = (d1 < d2) ? d2 : d1 ;
	  /* 0,0,0 and 1,1,1 */
	  like += (nn[0]+nn[7])*pStick[dmax] ;
	  /* 0,0,1 and 1,1,0 */
	  like += (nn[1]+nn[6])*(pStick[d1]+pSwitch[d2]-pSwitch[dmin]) ;
	  /* 0,1,0 and 1,0,1 */
	  like += (nn[2]+nn[5])*(pSwitch[d1]+pSwitch[d2]-pStick[dmin]) ;
	  /* 0,1,1 and 1,0,0 */
	  like += (nn[3]+nn[4])*(pSwitch[d1]+pStick[d2]-pSwitch[dmin]) ;
	}
      if (d1 < d2) { if (++d1 == d2) d2 = 0 ; } else { if (d2++ == d1) d1 = 0 ; }
    }
  
  free (pSwitch) ; free (pStick) ;

  return like ;
}

static double betaSearchLLDropOne (double beta)
{ return pbwtLogLikelihoodDropOne (info, alphaSearch, beta) ; }

static double alphaSearchLLDropOne (double alpha)
{ alphaSearch = alpha ;
  betaSearch = lineSearchPositive (betaSearch, 1.001, betaSearchLLDropOne) ;
  return pbwtLogLikelihoodDropOne (info, alphaSearch, betaSearch) ; 
}

/************************************************************/
/***** now a version using the column allele frequency ******/

static int pM ;

static Array buildRowInfoFreqDropOne (PBWT *p, int MAX) /* array of RowInfoDropOne */
/* record how many times we see each set of 3 consecutive values in y, 
   as a function of allele count */
{
  Array info = arrayCreate (p->M, RowInfoDropOne) ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  int i, j, k, n1 ;

  for (i = 0 ; i < p->N ; ++i)
    { for (j = 0 ; j < p->M ; ++j)
	{ if (!u->d[j] || !u->d[j+1]) continue ; /* ignore edge effects */
	  if (j == 0) 
	    k = (u->y[j] << 1) + u->y[j+1] ;
	  else if (j < p->M-1)
	    k = (u->y[j-1] << 2) + (u->y[j] << 1) + u->y[j+1] ;
	  else
	    k = (u->y[j-1] << 2) + (u->y[j] << 1) ;
	  n1 = u->M - u->c ;
	  arrayp(info,n1,RowInfoDropOne)->n[k] += 1 ; 
	  arrp(info,n1,RowInfoDropOne)->nTot += 1 ;
	}
      pbwtCursorForwardsReadAD (u, i) ;
    }

  pbwtCursorDestroy (u) ;

  pM = p->M ;

  return info ;
}

static double pbwtLLFreqDropOne (Array info, double alpha, double beta)
{
  int i ;
  double p00, p01, p10, p11, like = 0.0 ;
  RowInfoDropOne *inf = arrp(info,0,RowInfoDropOne) ;
  for (i = 0 ; i < arrayMax(info) ; ++i, ++inf) 
    if (inf->nTot)
      { double f = (0.5+i) / (double)(1+pM) ; /* frequency */
	p01 = -beta + alpha*log(f) ; p00 = log (1.0 - exp(p01)) ;
	p10 = -beta ; p11 = log (1.0 - exp(p10)) ;

	int *nn = inf->n ;
	/* 0,0,0 and 0,0,1 and 1,0,0 */ like += (nn[0]+nn[1]+nn[4])*p00 ;
	/* 0,1,0 */ like += nn[2]*(p01+p10-p00) ;
	/* 0,1,1 and 1,1,0 and 1,1,1 */ like += (nn[3]+nn[6]+nn[7])*p11 ;
	/* 1,0,1 */ like += nn[5]*(p10+p01-p11) ;
      }
  
  return like ;
}

static double betaSearchFreqDropOne (double beta)
{ return pbwtLLFreqDropOne (info, alphaSearch, beta) ; }

static double alphaSearchFreqDropOne (double alpha)
{ alphaSearch = alpha ;
  betaSearch = lineSearchPositive (betaSearch, 1.001, betaSearchFreqDropOne) ;
  return pbwtLLFreqDropOne (info, alphaSearch, betaSearch) ; 
}

/****************************************/

void pbwtFitAlphaBeta (PBWT *p, int model)
{
  double LL ;
  switch (model)
    {
    case 1:			/* alpha and beta drop one */
      info = buildRowInfoDropOne (p, 1000) ;
      alphaSearch = 0.0 ; 		/* first find beta-only model */
      betaSearch = lineSearchPositive (1.0, 1.001, betaSearchLLDropOne) ;
      LL = pbwtLogLikelihoodDropOne (info, alphaSearch, betaSearch) / p->N ;
      printf ("Fit beta %f  LL per site %f  per cell %f\n", betaSearch, LL, LL/p->M) ;
      alphaSearch = lineSearchPositive (0.01, 1.001, alphaSearchLLDropOne) ;
      LL = betaSearchLLDropOne (betaSearch) / p->N ;
      printf ("Fit alpha %f  beta %f  LL per site %f  per cell %f\n", 
	      alphaSearch, betaSearch, LL, LL/p->M) ;
      break ;
    case 2:			/* beta freq drop one */
      info = buildRowInfoFreqDropOne (p, 1000) ;
      alphaSearch = 1.0 ;
      betaSearch = lineSearchPositive (1.0, 1.001, betaSearchFreqDropOne) ;
      LL = pbwtLLFreqDropOne (info, alphaSearch, betaSearch) / p->N ;
      printf ("Fit beta %f  LL per site %f  per cell %f\n", betaSearch, LL, LL/p->M) ;
      alphaSearch = lineSearchPositive (1.0, 1.001, alphaSearchFreqDropOne) ;
      LL = betaSearchFreqDropOne (betaSearch) / p->N ;
      printf ("Fit alpha %f  beta %f  LL per site %f  per cell %f\n", 
	      alphaSearch, betaSearch, LL, LL/p->M) ;
      break ;
    }
  LL = -log(256.0)*arrayMax(p->yz) / p->N ;
  printf ("PBWT entropy per site %f  per cell %f\n", LL, LL/p->M) ;
  arrayDestroy (info) ;

  simpleEntropy (p) ;	/* print out simple entropy and some other stats */
}

/**************************************************************/
/******* Li and Stephens copying model ************************/

double copyLogLikelihoodDropOne (PBWT *p, double theta, double rho)
{
  int i, j, k ;
  PbwtCursor *u = pbwtCursorCreate (p, TRUE, TRUE) ;
  double **left = myalloc (p->M, double*) ;
  double *logLeftSum = mycalloc (p->M, double) ;
  for (i = 0 ; i < p->M ; ++i) 
    { left[i] = myalloc (p->M, double) ;
      for (j = 0 ; j < p->M ; ++j) left[i][j] = 1.0 / (p->M - 1.0) ;
      left[i][i] = 0.0 ;
    }
  uchar *x = myalloc (p->M, uchar) ;
  double rho1 = 1.0-rho, rhoM = rho/(p->M - 1.0), theta1 = 1.0-theta ;

  for (k = 0 ; k < p->N ; ++k)
    { for (j = 0 ; j < p->M ; ++j) x[u->a[j]] = u->y[j] ;
      for (i = 0 ; i < p->M ; ++i)
	{ double sum = 0.0 ;
	  for (j = 0 ; j < p->M ; ++j)
 	    { left[i][j] *= rho1 ;
	      left[i][j] += rhoM ;
	      left[i][j] *= (x[i] == x[j]) ? theta1 : theta ;
	      sum += left[i][j] ;
	    }
	  sum -= left[i][i] ; left[i][i] = 0.0 ;
	  logLeftSum[i] += log(sum) ;
	  for (j = 0 ; j < p->M ; ++j) if (j != i) left[i][j] /= sum ;
	}
      /*      if (isCheck) printf ("done site %d\n", k) ; */
      pbwtCursorForwardsRead (u) ;
    }
  pbwtCursorDestroy (u) ;

  double LL = 0 ;
  for (i = 0 ; i < p->M ; ++i) LL += logLeftSum[i] ;

  free (logLeftSum) ;
  for (i = 0 ; i < p->M ; ++i) free (left[i]) ; free (left) ;

  return LL ;
}

static PBWT *pSearch ;
double thetaSearch, rhoSearch ;

static double rhoSearchDropOne (double rho)
{ return copyLogLikelihoodDropOne (pSearch, thetaSearch, rho) ; }

static double thetaSearchDropOne (double theta)
{ thetaSearch = theta ;
  rhoSearch = lineSearchPositive (rhoSearch, 1.001, rhoSearchDropOne) ;
  return copyLogLikelihoodDropOne (pSearch, thetaSearch, rhoSearch) ; 
}

void pbwtLogLikelihoodCopyModel (PBWT *p, double theta, double rho)
{ double LL = copyLogLikelihoodDropOne (p, theta, rho) ;
  printf ("theta %f rho %f LL %f  per site %f  per cell %f\n", 
	  theta, rho, LL, LL/p->N, LL/(p->M*p->N)) ;
  pSearch = p ; 
  thetaSearch = theta ;
  rhoSearch = lineSearchPositive (rho, 1.01, rhoSearchDropOne) ;
  thetaSearch = lineSearchPositive (theta, 1.01, thetaSearchDropOne) ;
  LL = copyLogLikelihoodDropOne (pSearch, thetaSearch, rhoSearch) / p->N ;
  printf ("Fit theta %f  rho %f  LL per site %f  per cell %f\n", 
	  thetaSearch, rhoSearch, LL, LL/p->M) ;
}

/******** end of file ********/

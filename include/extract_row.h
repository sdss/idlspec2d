#ifndef __INCextract_row_h
#define __INCextract_row_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "extract_row.h"
  
#define TINY 1.0e-20

IDL_LONG extract_row (int argc, void *argv[]);
IDL_LONG extract_multi_row (int argc, void *argv[]);
void ProfileGauss(float *x, IDL_LONG ndat, float **y, float xcen, IDL_LONG xmin,
		IDL_LONG xmax, float sigma, IDL_LONG nCoeff);
void ProfileGaussStatic(float *x, IDL_LONG ndat, float **y, float xcen, 
                IDL_LONG xmin, IDL_LONG xmax, float sigma, IDL_LONG nCoeff);
void findXLimits(IDL_LONG *xmin, IDL_LONG *xmax, float *x, float *xcen, 
               IDL_LONG nx, IDL_LONG nTrace, float *sigma, float sigmal);
void fillProfile(float **y, float *x, float *xcen, IDL_LONG *xmin, 
             IDL_LONG *xmax, float *sigma, IDL_LONG nx, IDL_LONG nCoeff, 
             IDL_LONG nTrace, IDL_LONG proftype, float, float);
void fillPoly(float **y, float *x, IDL_LONG nx, IDL_LONG nPoly, 
         float x1, float x2);
void fillWhopping(float **y, float *x, IDL_LONG nx, IDL_LONG whoppingct,
         float *whoppingcen, float sigma);
void subtractProfile(float *y, IDL_LONG nx, IDL_LONG *xmin,
              IDL_LONG *xmax, IDL_LONG nTrace, IDL_LONG nCoeff,
              float **aprofile, IDL_LONG *ia, float *a);
void subtractPoly(float *y, IDL_LONG nx, IDL_LONG nPoly, float **apoly,
        IDL_LONG *ia, float *a);
void CheckRowFibers(float **abig, IDL_LONG *xmin, IDL_LONG *xmax, 
      IDL_LONG nTrace, IDL_LONG nCoeff, float *a, IDL_LONG *ia, float *invvar,
      float reject1, float reject2, float rejectflux, IDL_LONG *, IDL_LONG *);
void CheckMultiRowFibers(float **abig, IDL_LONG *xmin, IDL_LONG *xmax,
      IDL_LONG nTrace, IDL_LONG nCoeff, float *a, IDL_LONG *ia, float *invvar,
      IDL_LONG nsingle, IDL_LONG multirow);
void fillCovar(float *ysub, float *invvar, IDL_LONG nx, float **aprofile,
       float **apoly, IDL_LONG nTrace, IDL_LONG nCoeff, IDL_LONG nBand,
       IDL_LONG nPoly, float *beta, IDL_LONG *ia, float **covar, 
       IDL_LONG *xmin, IDL_LONG *xmax);
void cholslRow(float **a, IDL_LONG *ia, IDL_LONG nTrace, IDL_LONG nCoeff, 
         IDL_LONG nBand, IDL_LONG nPoly, float *p, float *b, float *x);
void cholslRowCovar(float **a, IDL_LONG *ia, IDL_LONG nTrace, IDL_LONG nCoeff, 
              IDL_LONG nBand, IDL_LONG nPoly, float *p);
int choldcRow(float **a, IDL_LONG *ia, IDL_LONG nTrace, IDL_LONG nCoeff, 
                 IDL_LONG nBand, IDL_LONG nPoly, float *p);
void chebyshevRow(float x, float *coeff, IDL_LONG nCoeff);


#endif /* __INCextract_row_h*/

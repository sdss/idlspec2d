
#ifndef __INCevilmath_h
#define __INCevilmath_h

void recenter_fweight
  (IDL_LONG    nx,
   float    *  imrow,
   float       radius,
   float    *  xcen,
   float    *  xerr);

int fixedGauss2(float x[], float y[], float ymod[], float invvar[], float *xcen,
        int ndat, int nTrace, int nPoly, float *sigma,
        float a[], int ia[], int ma, float **covar, float *chisq,
        int *xmin, int *xmax, float **abig);
float **matrix_nr(int nrow, int ncol);
void chebyshevFunc(float x, float *coeff, int nCoeff);

#endif /* __INCevilmath_h */


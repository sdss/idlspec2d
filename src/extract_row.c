#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "export.h"
#include "extract_row.h"
 
IDL_LONG extract_row
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   float     * x;
   float     * fimage;
   float     * invvar;
   float     * ymod;
   IDL_LONG    nTrace;
   IDL_LONG    nPoly;
   IDL_LONG    proftype;
   IDL_LONG    calcCovar;
   IDL_LONG    nCoeff;
   IDL_LONG    ma;
   float     * xcen;
   float     * sigma;
   IDL_LONG  * ia;
   float     * ans;
   float     * p;
   float    ** covar;
   float     * fscat;

   float     * whoppingcen;
   IDL_LONG    iy;
   IDL_LONG    whoppingct;
   IDL_LONG    squashprofile;
   IDL_LONG    retval = 1;
   IDL_LONG    bad;

   IDL_LONG    argct;
   IDL_LONG  * xmin;
   IDL_LONG  * xmax;
   float    ** aprofile;
   float    ** apoly;
   float       sigmal = 5.0; 	/* set limits of profile influence   */
   float       x2; 	/* Top Chebyshev x limit */
   float       x1; 	/* Lower Chebyshev x limit  */
   float       wsigma = 25.0; 	/* Whopping sigma width */
   IDL_LONG    wPoly; 	/* Whopping sigma width */

   IDL_LONG    i,j,k,l;
   IDL_LONG    length;
   IDL_LONG    mfit;
   IDL_LONG    coeff;
   IDL_LONG    tTrace;
   float     * ysub;
   float     * beta;

   /* Allocate pointers from IDL */

   argct  = 0;
   nx     = *((IDL_LONG *)argv[argct++]);
   x      = (float *)argv[argct++];
   fimage = (float *)argv[argct++];
   invvar = (float *)argv[argct++];
   ymod   = (float *)argv[argct++];

   nTrace = *((IDL_LONG *)argv[argct++]);
   nPoly  = *((IDL_LONG *)argv[argct++]);

   xcen   = (float *)argv[argct++];
   sigma  = (float *)argv[argct++];

   proftype = *((IDL_LONG *)argv[argct++]);
   calcCovar = *((IDL_LONG *)argv[argct++]);
   squashprofile = *((IDL_LONG *)argv[argct++]);
   whoppingcen = (float *)argv[argct++];
   whoppingct = *((IDL_LONG *)argv[argct++]);
   nCoeff = *((IDL_LONG *)argv[argct++]);
   ma     = *((IDL_LONG *)argv[argct++]);
   ans    = (float *)argv[argct++];
   ia     = (IDL_LONG *)argv[argct++];
   p      = (float *)argv[argct++];
   fscat  = (float *)argv[argct++];

   covar  = (float **)malloc(ma * sizeof(float *)); 
   for (iy=0; iy < ma; iy++) covar[iy] = (float *)argv[argct]+iy*ma;

/*   fprintf(stderr, "Going to fit_row now\n");  */

/*   for (i=0; i < nx; i++) 
      fprintf(stderr, "%d %f %f %f\n",(int) i,x[i],fimage[i],invvar[i]); */

   xmin = (IDL_LONG *)malloc(sizeof(IDL_LONG)*nTrace);
   xmax = (IDL_LONG *)malloc(sizeof(IDL_LONG)*nTrace);

   findXLimits(xmin, xmax, x, xcen, nx, nTrace, sigma, sigmal);

/*
//	ma = nCoeff*nTrace + nPoly
//	*/

   wPoly = nPoly + whoppingct;
   tTrace = nCoeff*nTrace;
   if (ma != tTrace + wPoly) {
      fprintf(stderr, "ma %d  does not equal nCoeff (%d) *nTrace (%d) + wPoly (%d)", (int) ma, (int) nCoeff, (int) nTrace, (int) wPoly);
      return -1;

   }	
   aprofile = (float **)malloc(tTrace * sizeof(float *));
   apoly = (float **)malloc(wPoly * sizeof(float *));

/*
//	Room for fiber profiles
//	*/

   for(i=0, k=0; i<nTrace; i++)
      for(j=0; j<nCoeff; j++, k++) {
         length = xmax[i] - xmin[i] + 1;
         if (length < 0) length = 0;
         aprofile[k] = (float *)malloc(length * sizeof(float));
      }

   fillProfile(aprofile, x, xcen, xmin, xmax, sigma, nx, nCoeff, 
          nTrace, proftype);

   if (squashprofile) {
/*     printf("Squashing Profile\n");  */
     for(i=0; i<nTrace; i++) {
       for(l=xmin[i],k=0;l<=xmax[i];l++,k++) {
	 aprofile[i*nCoeff][k] *= ans[i*nCoeff];
         for(j=1; j<nCoeff; j++) {
           aprofile[i*nCoeff][k] += aprofile[i*nCoeff+j][k]*ans[i*nCoeff+j];
	   }
      }
      for(j=1; j<nCoeff; j++) free(aprofile[i*nCoeff+j]);
     }

     for(i=0; i<nTrace; i++) {
       aprofile[i] = aprofile[i*nCoeff];
       ia[i] = 1;
/*       for(l=xmin[i],k=0;l<=xmax[i];l++,k++)  
//           printf("%f, ", aprofile[i][k]);
//	printf("\n%d",i+1); */
       }
     nCoeff = 1;
          
   }

/*
//	Room for polynomial profiles
//	*/

   for(i=0; i<wPoly; i++) {
      apoly[i] = (float *)malloc(nx * sizeof(float));
      ia[i+nCoeff*nTrace] = ia[i+tTrace];
      ans[i+nCoeff*nTrace] = ans[i+tTrace];
   }
    
   x2 = (float) (nx-1);
   x1 = 0.0;  
   if (nPoly > 0) fillPoly(apoly, x, nx, nPoly, x1, x2);

/*
//  Whopping profile has somewhere near 25 pixel sigma
//	*/
   fillWhopping(&apoly[nPoly], x, nx, whoppingct, whoppingcen, wsigma);

   CheckRowFibers(aprofile, xmin, xmax, nTrace, nCoeff, ans, ia, invvar);

/* Subtract out fixed variables first */

   ysub = (float *)malloc(nx * sizeof(float));

   for(i=0, mfit=0; i<ma; i++) if(ia[i]) mfit++;
   for(i=0; i<nx; i++) ysub[i] = fimage[i];


   if(mfit != ma && !squashprofile) {
/*      fprintf(stderr, "Subtracting fixed variables\n");  */
      subtractProfile(ysub, nx, xmin, xmax, nTrace, nCoeff, aprofile, ia, ans);
      subtractPoly(ysub, nx, wPoly, apoly, &ia[tTrace], &ans[tTrace]);
   } 

   beta = (float *)malloc(sizeof(float)*ma);

   for(i=0; i<ma; i++) 
      bzero((void *)covar[i] , sizeof(float)*ma);  

   for(i=0; i<ma; i++) beta[i] = 0.0;
   
   fillCovar(ysub, invvar, nx, aprofile, apoly, nTrace, nCoeff, wPoly, 
          beta, ia, covar, xmin, xmax);

   bad = choldcRow(covar, ia, nTrace, nCoeff, wPoly, p); 
/*   printf("choldc Custom2 done\n");	*/
   if (bad < 0) {
	for(j=0;j<ma;j++) {
	  p[j] = 0.0;
	  ans[j] = 0.0;
        }
   } else {

     cholslRow(covar, ia, nTrace, nCoeff, wPoly, p, beta, ans); 
/*     printf("cholsl done\n");	*/

     if (calcCovar > 0) {
        cholslRowCovar(covar, ia, nTrace, nCoeff, wPoly, p); 
/*     printf("cholsl Covar done\n");	*/
     }  
/*     else {
//      printf("Skipping Covariance Calculation\n");
//   }	*/
   }

   for(i=0;i<nx;i++) ymod[i] = 0.0;
    
   /* scattered light first  */

   for (j=nCoeff*nTrace,k=0;k<wPoly;j++,k++) {
      for(i=0; i < nx; i++) 
         ymod[i] += ans[j]*apoly[k][i];
    }
/*   printf("scattered light model done\n");  */

   /* Use 0th term to estimate scattered light contribution  */
   for (j=0;j<nTrace;j++) 
      for (i=xmin[j], k=0, fscat[j]=0.0; i <= xmax[j]; i++, k++) 
            fscat[j] += aprofile[j*nCoeff][k]*ymod[i];

   for (j=0,l=0;j<nTrace;j++) 
      for (coeff=0; coeff<nCoeff; coeff++,l++)
         for (i=xmin[j], k=0; i <= xmax[j]; i++, k++) 
            ymod[i] += ans[l]*aprofile[l][k];
/*   printf("model done\n");  */


   /* Free temporary memory */
   for(i=0,k=0; i<nTrace; i++) 
     for(j=0; j<nCoeff; j++,k++)
        free(aprofile[k]);
   for(i=0; i<wPoly; i++) free(apoly[i]);

   free(aprofile);
   free(apoly);
   free(ysub);
   free(beta);
   free(xmin);
   free(xmax);
   free(covar);
/*   printf("variables freed\n");   */

   return retval;
}

void ProfileGauss(float *x, IDL_LONG ndat, float **y, float xcen, IDL_LONG xmin,
		IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG nm,i,j,k,backup,place;
	float base;
	float diff, denom, frac;
	float sqbase, xfake;
	
	static float oldsigma=0.0;
	static float model[100][3][31];


	denom = 1.0/sqrt(6.2832 * sigma * sigma);

	if (sigma != oldsigma) {

	   printf("Filling Arrays\n");
	   oldsigma = sigma;
/*	   Fill static arrays	*/
	   for(nm=0;nm<100;nm++) {
 	     xfake = 15.0 + (float) nm/100.0;
	     for (k=0; k<=30; k++) {
	       for(j=0;j<3;j++) model[nm][j][k] = 0.0;
	       for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	     
	         diff = (xfake - (float)k + frac)/sigma;
	         base = exp(-diff*diff/2.0)*denom;

                 model[nm][0][k] += base;
                 sqbase = diff*diff*base;
	         model[nm][1][k] += sqbase;
	         model[nm][2][k] += diff*base;

                }
	     for (j=0;j<3;j++) model[nm][j][k] /= 5.0;
	   }
	}
      }

	frac = xcen - (int)xcen;
	nm = frac*100.0;
        backup = (int)xcen - xmin;
	for(j=0;j<nCoeff;j++) 
	  for (i=xmin,k=0,place=15-backup; i<=xmax; i++,k++,place++) 
	     { 
/*             printf("%d %d %d %d %d %f\n", 
		i,j,k,place,nm,frac);  */
             y[j][k] = model[nm][j][place];
	     }
	 
}

void ProfileAbs3(float *x, IDL_LONG ndat, float **y, float xcen, 
                IDL_LONG xmin, IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG nm,i,j,k,backup,place;
	float base;
	float diff, diffabs, denom, frac;
	float sqbase, xfake;
	
	static float oldsigma=0.0;
	static float model[100][3][31];

/*		Below is denominator for x^3	*/

	denom = 1.0/(2.88450 * 0.89298 * sigma);

	if (sigma != oldsigma) {

	   printf("Filling Arrays\n");
	   oldsigma = sigma;
/*	   Fill static arrays	*/
	   for(nm=0;nm<100;nm++) {
 	     xfake = 15.0 + (float) nm/100.0;
	     for (k=0; k<=30; k++) {
	       for(j=0;j<3;j++) model[nm][j][k] = 0.0;
	       for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	     
	         diff = (xfake - (float)k + frac)/sigma;
                 diffabs = fabs(diff);
	         base = exp(-diff*diff*diffabs/3.0)*denom;

                 model[nm][0][k] += base;
                 sqbase = diff*diff*diffabs*base;
	         model[nm][1][k] += sqbase;
	         model[nm][2][k] += diff*diffabs*base;

                }
	     for (j=0;j<3;j++) model[nm][j][k] /= 5.0;
	   }
	}
      }

	frac = xcen - (int)xcen;
	nm = frac*100.0;
        backup = (int)xcen - xmin;
	for(j=0;j<nCoeff;j++) 
	  for (i=xmin,k=0,place=15-backup; i<=xmax; i++,k++,place++) 
	     { 
/*             printf("%d %d %d %d %d %f\n", 
		i,j,k,place,nm,frac);	*/
             y[j][k] = model[nm][j][place];
	     }
	 
}

void ProfileAbs25(float *x, IDL_LONG ndat, float **y, float xcen, 
                IDL_LONG xmin, IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG nm,i,j,k,backup,place;
	float base;
	float diff, diffabs, denom2, frac;
	float denom3;
	float sqbase, xfake;
	
	static float oldsigma=0.0;
	static float model[100][3][31];

/*		Below is denominator for x^3	*/

	denom3 = 1.0/(2.88450 * 0.89298 * sigma);
	denom2 = 1.0/sqrt(6.2832 * sigma * sigma);

	if (sigma != oldsigma) {

	   printf("Filling Arrays\n");
	   oldsigma = sigma;
/*	   Fill static arrays	*/
	   for(nm=0;nm<100;nm++) {
 	     xfake = 15.0 + (float) nm/100.0;
	     for (k=0; k<=30; k++) {
	       for(j=0;j<3;j++) model[nm][j][k] = 0.0;
	       for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	     
	         diff = (xfake - (float)k + frac)/sigma;
                 diffabs = fabs(diff);
	         base = exp(-diff*diff*diffabs/3.0)*denom3;

/*	Add in ^3 term first  */
                 model[nm][0][k] += 1.5*base;
                 sqbase = diff*diff*diffabs*base;
	         model[nm][1][k] += 1.5*sqbase;
	         model[nm][2][k] += 1.5*diff*diffabs*base;

/*	Add in ^2 term next  */
	         base = exp(-diff*diff/2.0)*denom2;
                 model[nm][0][k] += 0.5*base;
                 sqbase = diff*diff*base;
	         model[nm][1][k] += 0.5*sqbase;
	         model[nm][2][k] += 0.5*diff*base;
                }
/*	Divide by 10 instead of 5 to get area 1.0  */
	     for (j=0;j<3;j++) model[nm][j][k] /= 10.0;
	   }
	}
      }

	frac = xcen - (int)xcen;
	nm = frac*100.0;
        backup = (int)xcen - xmin;
	for(j=0;j<nCoeff;j++) 
	  for (i=xmin,k=0,place=15-backup; i<=xmax; i++,k++,place++) 
	     { 
/*             printf("%d %d %d %d %d %f\n", 
		i,j,k,place,nm,frac);	*/
             y[j][k] = model[nm][j][place];
	     }
	 
}

void ProfileDoubleGauss(float *x, IDL_LONG ndat, float **y, float xcen, IDL_LONG xmin,
		IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG i,j,k;
	float base, mult;
	float diff, denom, frac;
        float sigma2, diff2, base2, denom2;
	float sqbase;

	sigma2 = 2.0*sigma;
	denom = 1.0/sqrt(6.2832 * sigma * sigma);
	denom2 = 1.0/sqrt(6.2832 * sigma2 * sigma2);


	for (i=xmin,k=0; i<=xmax; i++, k++) {
	   for (j=0;j<nCoeff;j++) y[j][k] = 0.0;
	  if(i >= 0 && i < ndat && nCoeff > 0) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = (xcen - x[i] + frac)/sigma;
	        diff2 = (xcen - x[i] + frac)/sigma2;
	        base = exp(-diff*diff/2.0)*denom;
	        base2 = exp(-diff2*diff2/2.0)*denom2;

                y[0][k] += base;
	        if(nCoeff >1) y[1][k] += base2;
	        if(nCoeff >2) {
                   sqbase = diff*diff*base;
	           y[2][k] += sqbase;

	           if(nCoeff >3) y[3][k] += diff*base;

	           for (j=4,mult=diff;j<nCoeff;j++,mult *= diff)
	              y[j][k] += mult*sqbase;
                   }
                }
	     for (j=0;j<nCoeff;j++) y[j][k] /= 5.0;
	   }
	}
}

void ProfileAbs3WideGauss(float *x, IDL_LONG ndat, float **y, float xcen, 
                IDL_LONG xmin, IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG nm,i,j,k,backup,place;
	float base, xfake;
	float diff, diffabs, denom, frac;
        float sigma2, diff2, base2, denom2;
	float sqbase;
	
	static float oldsigma=0.0;
	static float model[100][4][31];

/*		Below is denominator for x^3	*/

	denom = 1.0/(2.88450 * 0.89298 * sigma);
	sigma2 = 2.0*sigma;
	denom2 = 1.0/sqrt(6.2832 * sigma2 * sigma2);

	if (sigma != oldsigma) {

	   printf("Filling Arrays\n");
	   oldsigma = sigma;
/*	   Fill static arrays	*/
	   for(nm=0;nm<100;nm++) {
 	     xfake = 15.0 + (float) nm/100.0;
	     for (k=0; k<=30; k++) {
	       for(j=0;j<4;j++)model[nm][j][k] = 0.0;
	       for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	     
	         diff = (xfake - (float)k + frac)/sigma;
                 diffabs = fabs(diff);
	         base = exp(-diff*diff*diffabs/3.0)*denom;
	         diff2 = (xfake - (float)k + frac)/sigma2;
	         base2 = exp(-diff2*diff2/2.0)*denom2;


                 model[nm][0][k] += base;
	         model[nm][1][k] += base2;
                 sqbase = diff*diff*diffabs*base;
	         model[nm][2][k] += sqbase;
	         model[nm][3][k] += diff*diffabs*base;

                }
	     for (j=0;j<4;j++) model[nm][j][k] /= 5.0;
	   }
	}
      }

	frac = xcen - (int)xcen;
	nm = frac*100.0;
        backup = (int)xcen - xmin;
	for(j=0;j<nCoeff;j++) 
	  for (i=xmin,k=0,place=15-backup; i<=xmax; i++,k++,place++) 
	     {
/*             printf("%d %d %d %d %d %f\n", 
		i,j,k,place,nm,frac);	*/
             y[j][k] = model[nm][j][place];
	     }
	 
}


void ProfileAbs3HalfLorentz(float *x, IDL_LONG ndat, float **y, float xcen, 
                IDL_LONG xmin, IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG i,j,k;
	float base, mult;
	float diff, diffabs, denom, frac;
        float fwhm2, diff2, base2, denom2;
	float sqbase;
/*        float total; 	*/


/*		Below is denominator for x^3	*/

	denom = 1.0/(2.88450 * 0.89298 * sigma);
	fwhm2 = 2.0*sigma;
	denom2 = 3.1416/(2.0 * fwhm2);

	for (i=xmin,k=0; i<=xmax; i++, k++) {
	   for (j=0;j<nCoeff;j++) y[j][k] = 0.0;
	  if(i >= 0 && i < ndat && nCoeff > 0) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = (xcen - x[i] + frac)/sigma;
                diffabs = fabs(diff);
	        base = exp(-diff*diff*diffabs/3.0)*denom;
	        diff2 = 2.0*(xcen - x[i] + frac)/fwhm2;
	        base2 = denom2/(diff2*diff2 + 1.0);

                y[0][k] += base;
	        if(nCoeff >1 && diff2 < 0.0) y[1][k] += base2;
	        if(nCoeff >2 && diff2 >= 0.0) y[2][k] += base2;
	        if(nCoeff >3) {
                   sqbase = diff*diff*diffabs*base;
	           y[3][k] += sqbase;

	           if(nCoeff >4) y[4][k] += diff*base;

	           for (j=5,mult=diff;j<nCoeff;j++,mult *= diff)
	              y[j][k] += mult*sqbase;
                   }
                }
	     for (j=0;j<nCoeff;j++) y[j][k] /= 5.0;
	   }
	}

/*
//	This is just for debugging profiles
//
//	for (j=0;j<nCoeff;j++) { 
//	  for (i=xmin,k=0,total=0.0; i<=xmax; i++, k++) 
//	    if(i >= 0 && i < ndat && nCoeff > 0) {
//	      total += y[j][k];
//	    }
 //         fprintf(stderr,"%d %f ",(int)j,total);
//	}
 //         fprintf(stderr,"\n");	*/
	 
}


void ProfileGaussHalfLorentz(float *x, IDL_LONG ndat, float **y, float xcen, 
                IDL_LONG xmin, IDL_LONG xmax, float sigma, IDL_LONG nCoeff)
{ 
	IDL_LONG i,j,k;
	float base, mult;
	float diff, denom, frac;
        float fwhm2, diff2, base2, denom2;
	float sqbase;
/*        float total;	*/


/*		Below is denominator fro x^3	*/

	denom = 1.0/sqrt(6.2832 * sigma * sigma);
	fwhm2 = 2.0*sigma;
	denom2 = 5.0/(fwhm2 * 3.1416);

	for (i=xmin,k=0; i<=xmax; i++, k++) {
	   for (j=0;j<nCoeff;j++) y[j][k] = 0.0;
	  if(i >= 0 && i < ndat && nCoeff > 0) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = (xcen - x[i] + frac)/sigma;
	        base = exp(-diff*diff/2.0)*denom;
	        diff2 = 2.0*(xcen - x[i] + frac)/fwhm2;
	        base2 = denom2/(diff2*diff2 + 1.0);

                y[0][k] += base;
	        if(nCoeff >1 && diff2 < 0.0) y[1][k] += base2;
	        if(nCoeff >2 && diff2 >= 0.0) y[2][k] += base2;
	        if(nCoeff >3) {
                   sqbase = diff*diff*base;
	           y[3][k] += sqbase;

	           if(nCoeff >4) y[4][k] += diff*base;

	           for (j=5,mult=diff;j<nCoeff;j++,mult *= diff)
	              y[j][k] += mult*sqbase;
                   }
                }
	     for (j=0;j<nCoeff;j++) y[j][k] /= 5.0;
	   }
	}

/*
//	This is just for debugging profiles
//
//	for (j=0;j<nCoeff;j++) { 
//	  for (i=xmin,k=0,total=0.0; i<=xmax; i++, k++) 
//	    if(i >= 0 && i < ndat && nCoeff > 0) {
//	      total += y[j][k];
//	    }
 //         fprintf(stderr,"%d %f ",(int)j,total);
//	}
 //         fprintf(stderr,"\n");	*/
	 
}


void findXLimits(IDL_LONG *xmin, IDL_LONG *xmax, float *x, float *xcen, 
               IDL_LONG nx, IDL_LONG nTrace, float *sigma, float sigmal)
{
/*
//        Find max and min of profile influence
//	*/
   IDL_LONG i, place, place2;
   float diff;

   place = 0;

   for(i=0; i< nTrace; i++) {
      diff = sigmal*sigma[i];
      while(x[place] < xcen[i] - diff && place < nx) place++;
      if(place >= nx) place--;

      place2 = place;
      while(x[place2] <= xcen[i] + diff && place2 < nx) place2++;
      if(place2 >= nx) place2--;

      xmin[i] = place;
      xmax[i] = place2 - 1; 

/*      if(xmax[i] <= xmin[i]) 
         fprintf(stderr," Fiber # %d has no influence %f \n", (int) i, diff); */

      place = xmin[i];
   }     

}

void subtractProfile(float *y, IDL_LONG nx, IDL_LONG *xmin, 
              IDL_LONG *xmax, IDL_LONG nTrace, IDL_LONG nCoeff,
              float **aprofile, IDL_LONG *ia, float *a)
{
	IDL_LONG i,j,k,l,coeff;

        for (i=0,l=0;i<nTrace;i++) 
           for (coeff=0; coeff<nCoeff; coeff++,l++) 
              if (!ia[l])  
                 for (j=xmin[i],k=0;j<=xmax[i];j++,k++)  
                    y[j] -= a[l] * aprofile[l][k];
               
}

void subtractPoly(float *y, IDL_LONG nx, IDL_LONG nPoly, float **apoly, 
        IDL_LONG *ia, float *a)
{
	IDL_LONG i,j;
        
        for (i=0;i<nPoly;i++) 
           if (!ia[i]) 
              for (j=0;j<nx;j++) 
                 y[i] -= a[i] * apoly[i][j];
}
           
void CheckRowFibers(float **abig, IDL_LONG *xmin, IDL_LONG *xmax, 
      IDL_LONG nTrace, IDL_LONG nCoeff, float *a, IDL_LONG *ia, float *invvar) 
{

	IDL_LONG i,j,k,l,m;
	float total;
   
	for (i=0;i<nTrace;i++) {
	   total = 0.0;
	   for (k=xmin[i],m=0; k<=xmax[i]; k++,m++)
	      if (invvar[k] > 0.0) total += abig[i*nCoeff][m];
	   if (total < 0.8) {
              for(j=nCoeff-1,l=j+i*nCoeff;j>=1;j--,l--) 
	         if (ia[l]) {
	            ia[l] = 0;
	            a[l] = 0.0;
/*          fprintf(stderr,"Fiber %d, dropped term %d, total: %f\n", 
                    (int) i, (int) j,total);   */
                    }
                 }
	   if (total < 0.4) {
              j = 0;
              l=i*nCoeff;
	      if (ia[l]) {
	         ia[l] = 0;
	         a[l] = 0.0;
/*          fprintf(stderr,"Fiber %d, dropped term %d, total: %f\n", 
                   (int) i, (int) j,total);   */
                    }
               }
	}   
}

void CheckMultiRowFibers(float **abig, IDL_LONG *xmin, IDL_LONG *xmax, 
      IDL_LONG nTrace, IDL_LONG nCoeff, float *a, IDL_LONG *ia, float *invvar, 
      IDL_LONG nsingle, IDL_LONG multirow) 
{

	IDL_LONG i,j,k,l,m;
	IDL_LONG row;
	float total;
   
	for (i=0;i<nTrace;i++) {
	   total = 0.0;
           for (row = 0; row < multirow; row++)
	   for (k=xmin[i+row*nTrace],m=0; k<=xmax[i+row*nTrace]; k++,m++)
	      if (invvar[k+row*nsingle] > 0.0) total += abig[(i+row*nTrace)*nCoeff][m];
	         if (total < 0.8) {
                    for(j=nCoeff-1,l=j+i*nCoeff;j>=1;j--,l--) 
	               if (ia[l]) {
	                  ia[l] = 0;
	                   a[l] = 0.0;
/*          fprintf(stderr,"Fiber %d, dropped term %d, total: %f\n", 
                    (int) i, (int) j,total);   */
                    }
                 }
	      if (total < 0.4) {
                 j = 0;
                 l=i*nCoeff;
	            if (ia[l]) {
	               ia[l] = 0;
	                a[l] = 0.0;
/*          fprintf(stderr,"Fiber %d, dropped term %d, total: %f\n", 
                   (int) i, (int) j,total);   */
                    }
               }
	}   
}

void fillProfile(float **y, float *x, float *xcen, IDL_LONG *xmin, 
             IDL_LONG *xmax, float *sigma, IDL_LONG nx, IDL_LONG nCoeff, 
             IDL_LONG nTrace, IDL_LONG proftype)
{
	int i, j, length;

	for (i=0,j=0; i<nTrace; i++,j+=nCoeff) {
	    length = xmax[i] - xmin[i] + 1;
	    if (length > 0 && sigma[i] > 0.0) {
	      if (proftype == 1) 
                 ProfileGauss(x, nx, &y[j], xcen[i], xmin[i], xmax[i], 
                      sigma[i], nCoeff);
	      else if (proftype == 2)  {
                 ProfileAbs3(x, nx, &y[j], xcen[i], xmin[i], xmax[i], 
                      sigma[i], nCoeff);
                 }
	      else if (proftype == 3)  {
                 ProfileAbs25(x, nx, &y[j], xcen[i], xmin[i],xmax[i],
                      sigma[i], nCoeff);
                 }
	      else if (proftype == 4)  {
                 ProfileDoubleGauss(x, nx, &y[j], xcen[i], xmin[i],xmax[i],
                      sigma[i], nCoeff);
                 }
	      else if (proftype == 5)  {
                 ProfileAbs3WideGauss(x, nx, &y[j], xcen[i], xmin[i],xmax[i],
                      sigma[i], nCoeff);
                 }
	      else if (proftype == 6)  {
                 ProfileAbs3HalfLorentz(x, nx, &y[j], xcen[i], xmin[i],xmax[i],
                      sigma[i], nCoeff);
                 }
	      else if (proftype == 7)  {
                 ProfileGaussHalfLorentz(x, nx, &y[j], xcen[i], xmin[i],xmax[i],
                      sigma[i], nCoeff);
                 }
              else {
	         fprintf(stderr,"Using Gaussian");
                 ProfileGauss(x, nx, &y[j], xcen[i], xmin[i], xmax[i], 
                      sigma[i], nCoeff);
              }
            }
	}
 
}

void fillPoly(float **y, float *x, IDL_LONG nx, IDL_LONG nPoly, 
         float x1, float x2) 
{
	int i,j;
        float norm;

	float *atemp = (float *)malloc(nPoly * sizeof(float));   

        for (i=0;i<nx;i++) {
           norm = (2.0 * x[i] - (x1 + x2))/(x2 - x1);

/*	These next two lines make a step function to account for
	uncertainties in the amplifier's gains matching 

	   if (norm < 0.0) y[0][i] = 0.0; 
	   else y[0][i] = 1.0; 

	The rest of the params are just nPoly-1 chebyshev coefficients 

	   chebyshevRow(norm, atemp, nPoly-1);
	   for(j=1; j<nPoly; j++) {
	     y[j][i] = atemp[j-1];
           }   

	     I don't think the step function is helping, it seems
	     to be giving wrong answers when traces fall off one side */


	   chebyshevRow(norm, atemp, nPoly);
	   for(j=0; j<nPoly; j++) {
	     y[j][i] = atemp[j];
           }


        }

	free(atemp);
}

void fillWhopping(float **y, float *x, IDL_LONG nx, IDL_LONG whoppingct, 
         float *whoppingcen, float sigma) 
{
	int i,j;
        float diff;

        for (j=0;j<whoppingct;j++) {
          for (i=0;i<nx;i++) {
	     diff = fabs(whoppingcen[j] - x[i])/sigma;
	     y[j][i] = 0.5*exp(-diff)/sigma;
             }
          }

}

void fillCovar(float *ysub, float *invvar, IDL_LONG nx, float **aprofile, 
       float **apoly, IDL_LONG nTrace, IDL_LONG nCoeff, IDL_LONG nPoly, 
       float *beta, IDL_LONG *ia, float **covar, IDL_LONG *xmin, IDL_LONG *xmax)
{ 

   IDL_LONG i,j,k,l,m,n;
   IDL_LONG tTrace = nCoeff*nTrace;
   IDL_LONG coeff;
   IDL_LONG mStop,jStop;

/* Fill Profile part of covar first  */

  for (l=0,j=0,jStop = nCoeff,mStop=2*nCoeff;
             l<nTrace;l++,mStop+=nCoeff,jStop+=nCoeff) {
    if(mStop > tTrace) mStop = tTrace;
    for (coeff=0; coeff < nCoeff; coeff++, j++)
       if (ia[j]) {
         for (i=xmin[l],k=0;i<=xmax[l];i++,k++) 
            beta[j] += ysub[i] * aprofile[j][k] * invvar[i];
         
	 for (m=j; m<jStop; m++) 
            if (ia[m]) 
               for (i=xmin[l],k=0;i<=xmax[l];i++,k++) 
                  covar[j][m] += aprofile[m][k] * aprofile[j][k] * invvar[i];
         

/*	Overlap with polynomial  */

         for (k=0,m=tTrace; k<nPoly; k++, m++)
            if (ia[m]) 
               for (i=xmin[l],n=0;i<=xmax[l];i++,n++) 
                  covar[j][m] += apoly[k][i] * aprofile[j][n] * invvar[i];
         

/*	Overlap with next fiber */

	 for (m=jStop; m<mStop; m++) 
            if (ia[m]) 
               for (i=xmin[l+1],k=0;i<=xmax[l];i++,k++) 
                 covar[j][m] += aprofile[m][k] * 
                    aprofile[j][i-xmin[l]] * invvar[i];
      }
      
   }

   for (l=0, j=tTrace; l <nPoly; l++, j++)
      if (ia[j]) {
         for (i=0; i<nx; i++) 
            beta[j] += ysub[i] * apoly[l][i] * invvar[i];

         for (k=l,m=j; k<nPoly; k++, m++)
            if (ia[m]) 
               for (i=0; i<nx; i++) 
                  covar[j][m] += apoly[k][i] * apoly[l][i] * invvar[i];
      }
	
/*   printf("nPoly done\n");	*/
}

/* cholslRow replaces lower triangle of a with sqrt(covar)  
	And we only loop over enough to work on nCoeff parameters for
	adjacent fibers   
   a is modified and cannot be used again for subsequent calls */
void cholslRow(float **a, IDL_LONG *ia, IDL_LONG nTrace, IDL_LONG nCoeff, 
         IDL_LONG nPoly, float *p, float *b, float *x) 
{
	IDL_LONG i,j,k;
        IDL_LONG kStart;
	IDL_LONG coeff;
	IDL_LONG tTrace=nTrace*nCoeff;
	float sum;

	for(j=0,i=0;j<nTrace;j++) {
          kStart = i - 2*nCoeff;
          if (kStart < -1) kStart = -1;
          for(coeff=0;coeff<nCoeff;coeff++,i++)
             if (ia[i]) {
	        sum = b[i];
	   
	      for (sum=b[i],k=i-1; k > kStart; k--) 
                 if (ia[k]) sum -= a[i][k]*x[k];
	      x[i] = sum/p[i];
	   }
        }
	for(i=tTrace;i<tTrace+nPoly;i++) 
           if (ia[i]) {
	      for (sum=b[i],k=i-1;k>=0;k--) 
              if (ia[k])  sum -= a[i][k]*x[k];
	      x[i] = sum/p[i];
	}


	for(i=tTrace + nPoly - 1;i>=tTrace;i--) 
           if (ia[i]) {
	      for (sum=x[i],k=i+1;k<tTrace+nPoly;k++) 
                 if (ia[k])  sum -= a[k][i]*x[k];
	   x[i] = sum/p[i];
	}

	for (j=nTrace-1,i=tTrace-1,kStart=tTrace;j>=0;j--) {
           if (j < nTrace-2) kStart -= nCoeff;
           for(coeff=0;coeff<nCoeff;coeff++,i--)
           if (ia[i]) {
	      for (sum=x[i],k=i+1; k<kStart; k++) 
                 if (ia[k])  sum -= a[k][i]*x[k];
	      for (k=tTrace;k<tTrace + nPoly;k++) 
                 if (ia[k]) sum -= a[k][i]*x[k];
	      x[i] = sum/p[i];
	   }
        }
}

void cholslRowCovar(float **a, IDL_LONG *ia, IDL_LONG nTrace, IDL_LONG nCoeff, 
              IDL_LONG nPoly, float *p)
{
	IDL_LONG i,k;
        IDL_LONG kStart = 2*nCoeff;
	IDL_LONG tTrace=nTrace*nCoeff;
	float sum;
	float *x;
	IDL_LONG pl;
	
	for(pl=0;pl < tTrace+nPoly; pl++)  
           if (ia[pl]) {
	      x = a[pl];

	      for(i=pl,sum=1.0;i<tTrace;i++,sum=0.0) 
                 if (ia[i]) {
	         for (k=i-1;k>=pl && k>=i-kStart;k--) 
                    if (ia[k]) sum-=a[i][k]*x[k];
	         x[i] = sum/p[i];
	      }
	      for(;i<tTrace+nPoly;i++,sum=0.0) 
                 if (ia[i]) {
                    for (k=i-1;k>=pl;k--) 
                       if (ia[k]) sum -= a[i][k]*x[k];
	         x[i] = sum/p[i];
	      }

	      for (i=tTrace+nPoly-1; i>=tTrace && i >= pl;i--) 
                 if (ia[i]) {
                    sum = x[i];
	            for (k=i+1;k<tTrace+nPoly;k++) 
                       if (ia[k]) sum -= a[k][i]*x[k];
	            x[i] = sum/p[i];
	      }

	      for(i=tTrace-1;i>=0 && i >= pl;i--) 
                 if (ia[i]) {
                    sum = x[i];
	            for (k=i+1;k<tTrace && k<i+kStart;k++) 
                       if (ia[k]) sum -= a[k][i]*x[k];
	            for (k=tTrace;k<tTrace + nPoly;k++) 
                       if (ia[k]) sum -= a[k][i]*x[k];
	            x[i] = sum/p[i];
	      }
	  }

}

/*	Row requires nCoeff  parameters per Trace in upper triangle of covar */           
int choldcRow(float **a, IDL_LONG *ia, IDL_LONG nTrace, IDL_LONG nCoeff, 
                 IDL_LONG nPoly, float *p)
{
   IDL_LONG i,j,k,l;
   IDL_LONG coeff;
   IDL_LONG tTrace = nCoeff*nTrace;
   IDL_LONG kStart;
   IDL_LONG jStop;
   float sum;

   for (l=0,i=0,jStop=2*nCoeff;l<nTrace;l++,jStop+=nCoeff) {
      if (jStop > tTrace) jStop = tTrace;
      kStart = i - 2*nCoeff;
      if(kStart < -1) kStart = -1;
      for (coeff=0; coeff < nCoeff; coeff++, i++)
         if (ia[i]) {
            for (j=i; j<jStop; j++) 
               if (ia[j]) {
                  for (sum = a[i][j], k=i-1; k > kStart; k--) 
	             if(ia[k]) sum -= a[i][k]*a[j][k];

                     if (i==j) {
                        if (sum <= 0.0) {
                           fprintf(stderr,"choldc failed %d %f\n",(int) i,sum);
                           sum = a[i][j];
                           fprintf(stderr,"%f\n",sum);
                           for (k=i-1; k > kStart; k--) {
                              if(ia[k]) sum -= a[i][k]*a[j][k];
                           fprintf(stderr,"%f %f %f\n",sum, a[i][k], a[j][k]);
                        }
                     return -i;
                     }
                  p[i] = sqrt(sum);
                  }
               else a[j][i] = sum/p[i];
            }

	    for(j=tTrace; j<tTrace +nPoly; j++) 
	       if(ia[j]) {
                  for(sum = a[i][j], k=i-1; k > kStart; k--) 
	             if(ia[k]) sum -= a[i][k]*a[j][k];
              
                  a[j][i] = sum/p[i];
	       }
            }
        }

	for(i=tTrace;i<tTrace+nPoly;i++) 
	   if(ia[i]) {
	   for(j=i;j<tTrace+nPoly;j++) 
	      if(ia[j]) {
	         for (sum=a[i][j],k=i-1;k>=0;k--) 
	            if(ia[k]) sum -= a[i][k]*a[j][k];
	         if(i==j) {
	            if (sum <= 0.0)  {
                      fprintf(stderr,"choldc failed %d\n", (int) i);
                      return -i;
                    }
	            p[i] = sqrt(sum);
	      } else a[j][i] = sum/p[i];
           }
        }
	return 0;
}

void chebyshevRow(float x, float *coeff, IDL_LONG nCoeff)
{
        /* Return Terms of Chebyshev Polynomial */

        IDL_LONG i;
        float twox;

        if(nCoeff > 0) coeff[0] = 1.0;
        if(nCoeff > 1) coeff[1] = x;

        twox = 2.0*x;
        for(i=2;i<nCoeff;i++)
           coeff[i] = twox * coeff[i-1] - coeff[i-2];

      return;
}


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "evilmath.h"
  

void fillAnswers(float *a, int *ia, float **covar, int nTrace, int nPoly, 
                   float *afit, float *error, float *sigma, float *sigmaerr)
{
	int i,j,change;
	float sum;
	float covar00, covar11, covar01;
	float r;
	float *sigmanew = (float *)malloc(nTrace * sizeof(float *));
	
	sum = 0.0;
	j = -1;
	change = 0;

	for(i=0;i<nTrace;i++) {
	  sigmanew[i] = sigma[i];
	  sigmaerr[i] = 0.0;

	  afit[i] = a[i*2] + a[i*2+1];
	  covar00 = 0.0;
	  covar01 = 0.0;
	  covar11 = 0.0;

	  if(ia[i*2]) {
	    j++;
	    covar00 = covar[j][j];
	    if(ia[i*2+1]) covar01 = covar[j][j+1];
	  }
	  if(ia[i*2+1]) {
	    j++;
	    covar11 = covar[j][j];
	  }
	 
	  error[i] = sqrt(covar00 + covar11 + 2.0*covar01);

	  if (a[i*2] > 0.0) {
	     r = a[i*2+1]/a[i*2];
	     sigmanew[i] = sigma[i] * (1.0 + r);
	     sigmaerr[i] = sigma[i]/a[i*2] *
	        sqrt(covar00 + covar11 * r *r - 2.0 * covar01);
	     if (sigmanew[i] > 0.0 && sigmaerr[i] > 0.0) 
	        if (sigmaerr[i] / sigmanew[i] < 0.05)  {
                   change++;	
//	           printf("%d Changing width of fiber %d %f -> %f \n", 
//	           change, i, sigma[i], sigmanew[i]);
	           sigma[i] = sigmanew[i];
	        }
	     }
        }

/*	for (i=nTrace;i<nTrace+nPoly;i++) {
	   afit[i+nTrace] = a[i+nTrace];
	   error[i+nTrace] = 0.0;
	   if (ia[i+nTrace]) {
	      j++; 
	      error[i+nTrace] = sqrt(covar[j][j]);
	   }
	}  */

	free(sigmanew);
}

int fit_row(float *x, float *y, float *invvar, float *ymod, int nx,
              float *xcen, float *sigma, int nTrace, int nPoly, 
 	      float *fextract, float *ferror, float *fscattered, float *fwidth, 
              int maxIter, int refit, float highrej, float lowrej, float boxap)
{

	int i,j,k;
	int totalp;
	float *a;
        float **covar, **abig;
	float chisq, diffchi, oldchi;
	int rejected;
        int *xmin, *xmax;
        int *ia;
	int retval;
        float *sigmaerr = (float *)malloc(nTrace * sizeof(float *));
	float upperl, lowerl;
	float diff;

        totalp = 2*nTrace + nPoly;
        fprintf(stderr, "Inside fit_row %d %f %d %d %d %d %f %f %f\n", 
              nx, sigma[0], nTrace, 
              nPoly, maxIter, refit, highrej, lowrej, boxap);

 	a = (float *)malloc(totalp * sizeof(float)); 
 	ia = (int *)malloc(totalp * sizeof(int)); 

        xmin = (int *)malloc(nTrace * sizeof(int));
        xmax = (int *)malloc(nTrace * sizeof(int));
//        abig = (float **)malloc(totalp * sizeof(float *));
        covar = matrix_nr(totalp,totalp);
        abig = matrix_nr(totalp,20);


	for(i=0;i<nTrace;i++) fwidth[i] = sigma[i];
	for(i=0;i<totalp;i++) ia[i] = 1;

//	for(i=0;i<nx;i++) 
//          fprintf(stderr, "%f %f %f %f\n", x[i], y[i], ymod[i], invvar[i]);


        diffchi = 10.0;
        oldchi = 1.0e10;
	rejected = 10;
        k=0;
        while((diffchi > 0.5 || rejected > 0) && k < maxIter) {

          retval = fixedGauss2(x, y, ymod, invvar, xcen, nx, nTrace,
               nPoly, fwidth, a, ia, totalp, covar, &chisq, xmin, xmax, abig);
  
           fprintf(stderr, "Hello %d\n",maxIter);
          fillAnswers(a,ia,covar,nTrace,nPoly,fextract,ferror,fwidth,sigmaerr); 

	  for(j=0;j<nx;j++) fscattered[j] = 0.0;
          for(i=2*nTrace;i<2*nTrace+nPoly;i++) {
            printf("Background %d %f \n", i - 2*nTrace, a[i]);
	    for(j=0;j<nx;j++) fscattered[j] += abig[i][j]*a[i];
          }


        upperl = highrej*sqrt(chisq/(nx-totalp));
        lowerl = -lowrej*sqrt(chisq/(nx-totalp));
	rejected = 0;
        for(i=0;i<nx;i++) {
          diff = (y[i]-ymod[i]) * sqrt(invvar[i]);
          if (diff > upperl || diff < lowerl) {
            printf("Rejected %d %5.3f %5.3f %8.3f %8.3f %8.3f \n",
		i,upperl, lowerl,diff,y[i],ymod[i]);
            invvar[i] = 0.0;
	    rejected++;
          }
        }

        diffchi = oldchi - chisq;
        oldchi = chisq;
	k++;
      }

	free(a);
	free(ia);
	free(xmin);	
	free(xmax);	
	free(sigmaerr);	
        for(i=0;i<totalp;i++) free(abig[i]);
        for(i=0;i<totalp;i++) free(covar[i]);
	free(abig);	
	free(covar);	

	return retval;
}

IDL_LONG extract_profile
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float    ** invvar;
   float    ** ymod;
   IDL_LONG    nTrace;
   IDL_LONG    nRowExtract;
   IDL_LONG    nPoly;
   float    ** xcen;
   IDL_LONG  * ycen;
   float    ** sigma;
   float    ** fextract;
   float    ** ferror; 
   float    ** fscattered;
   float    ** fwidth;

   long        iy;
   long        ix;
   float    *  xdummy;
   IDL_LONG    retval = 1;

   IDL_LONG    maxIter;
   IDL_LONG    refit;
   float       highrej;
   float       lowrej;
   float       boxap;
   int         rowExtract;
   int         j;

   /* Allocate pointers from IDL */

   nx = *((IDL_LONG *)argv[0]);
   ny = *((IDL_LONG *)argv[1]);
   fimage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) fimage[iy] = (float *)argv[2] + iy*nx;

   invvar = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) invvar[iy] = (float *)argv[3] + iy*nx;

   ymod = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) ymod[iy] = (float *)argv[4] + iy*nx;

   nTrace = *((IDL_LONG *)argv[5]);
   nRowExtract = *((IDL_LONG *)argv[6]);

   xcen = (float **)malloc(nRowExtract * sizeof(float *)); 
   for (iy=0; iy < nRowExtract; iy++) xcen[iy] = (float *)argv[7] + iy*nTrace;

   ycen = (IDL_LONG *)argv[8]; 

   sigma = (float **)malloc(nRowExtract * sizeof(float *)); 
   for (iy=0; iy < nRowExtract; iy++) sigma[iy] = (float *)argv[9] + iy*nTrace;

   /* fextract has 5 rows now */
   fextract = (float **)malloc(nRowExtract * sizeof(float *)); 
   for (iy=0; iy < nRowExtract; iy++) 
      fextract[iy] = (float *)argv[10] + iy*nTrace;

   ferror = (float **)malloc(nRowExtract * sizeof(float *)); 
   for (iy=0; iy < nRowExtract; iy++) 
      ferror[iy] = (float *)argv[11] + iy*nTrace;

   fscattered = (float **)malloc(nRowExtract * sizeof(float *)); 
   for (iy=0; iy < nRowExtract; iy++) 
      fscattered[iy] = (float *)argv[12] + iy*nx;  
/*  Change fscattered back to *nTrace when done!!!! */

   fwidth= (float **)malloc(nRowExtract * sizeof(float *)); 
   for (iy=0; iy < nRowExtract; iy++) 
      fwidth[iy] = (float *)argv[13] + iy*nTrace;

   nPoly   = *((IDL_LONG *)argv[14]);
   maxIter = *((IDL_LONG *)argv[15]);
   refit   = *((IDL_LONG *)argv[16]);
   highrej = *((float *)argv[17]);
   lowrej  = *((float *)argv[18]);
   boxap   = *((float *)argv[19]);

   /* make dummy x array */
     xdummy = (float *)malloc(nx * sizeof(float)); 
     for (ix=0; ix < nx; ix++) xdummy[ix] = ix;


   /* Loop through each row */
   for (iy=0; iy< nRowExtract; iy++) {

     rowExtract = ycen[iy];
      fprintf(stderr, "%d %d\n",(int)iy, rowExtract);

   /* Initialize Extraction output */
     for(j=0; j<nTrace; j++) {
	fextract[iy][j] = 0.0;
	ferror[iy][j] = -1.0;
	fscattered[iy][j] = 0.0;
	fwidth[iy][j] = 0.0;
     }	
     
     if(rowExtract >= 0 && rowExtract < ny) { 
     
      /* Fit profiles at all trace centers, include some background
        if necessary */

        fprintf(stderr, "Going to fit_row\n");

	retval = fit_row(xdummy, fimage[rowExtract],invvar[rowExtract], 
                  ymod[rowExtract], (int) nx,
		  xcen[iy], sigma[iy], (int) nTrace, (int) nPoly, 
		  fextract[iy], ferror[iy], fscattered[iy], fwidth[iy], 
                  (int) maxIter, (int) refit, highrej, lowrej, boxap);
	}
   }

   /* Free temporary memory */
   free(xdummy);
   free(fimage);
   free(invvar);
   free(ymod);
   free(xcen);
   free(sigma);
   free(fextract);
   free(ferror);
   free(fscattered);
   free(fwidth);

   return retval;
}


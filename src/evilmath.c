
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "export.h"
#include "evilmath.h"

#define TINY 1.0e-20
 
/* Recenter one position in a row of data using a flux-weighted center */
void recenter_fweight
  (IDL_LONG    nx,
   float    *  imrow,
   float    *  invvar,
   float       radius,
   float       xinit,
   float    *  xcen,
   float    *  xerr)
{
   IDL_LONG    ii;
   IDL_LONG    ix1;
   IDL_LONG    ix2;
   IDL_LONG    npix;
   int         qbad;
   float       x1;
   float       x2;
   float       sumxw;
   float       sumf;
   float       sumsxsx;
   float       sumw;
   float       xdiff;
   float    *  convol;

   /* Only recenter if the guess value is within the bounds of the image */
   if (xinit > 0.0 && xinit < nx-1) {

      /* Determine which pixel numbers over which to sum */
      x1 = xinit - radius + 0.5;
      x2 = xinit + radius + 0.5;
      ix1 = floor(x1);
      ix2 = floor(x2);

      /* If either end of the convol is out of bounds in the image, then
         shrink both sides of the summing convol until it is in bounds.
       */
      if (x1 < 0.0) {
         x2 += x1;
         x1 = 0.0;
      }
      if (x2 > nx) {
         x1 += x2 - nx;
         x2 = nx;
         ix2 = nx - 1;
      }
      npix = ix2 - ix1 + 1;

      convol = (float *) malloc(sizeof(float) * npix);
      sumw = 0.0;
      sumxw = 0.0;
      sumf = 0.0;
      qbad = 0;
      for (ii=0; ii < npix; ii++) {
         /* Determine the weight of a boxcar window function (convol) for this
          * pixel.  Note that the values of "convol" will sum to 2*radius
          * unless the edge of the image is reached.
          */
         if (ii == 0) {
            convol[ii] = 1.0 + ix1 - x1;
         } else if (ii == npix-1) {
            convol[ii] = x2 - ix2;
         } else {
            convol[ii] = 1.0;
         }
         sumw = sumw + convol[ii] * imrow[ix1+ii];

         xdiff = ix1 + ii - xinit;
         sumxw += convol[ii] * imrow[ix1+ii] * xdiff;

         if (invvar[ix1+ii] <= 0.0) qbad = 1;
      }

      if (sumw > 0.0 && qbad == 0) {
         *xcen = sumxw / sumw + xinit;

         /* Compute the error in the flux-weighted center */
         sumsxsx = 0.0;
         for (ii=0; ii < npix; ii++) {
            xdiff = ix1 + ii - *xcen;
            sumsxsx += xdiff * xdiff * convol[ii] * convol[ii]
             / (sumw * sumw * invvar[ix1+ii]);
         }
         *xerr = sqrt(sumsxsx);
      } else {
         *xcen = xinit;
         *xerr = 999.0;
      }

      free(convol);
   }
}

float **matrix_nr(int nrow, int ncol)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        int i;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc(nrow*sizeof(float*));

        /* allocate rows and set pointers to them */
        for(i=0;i<nrow;i++) 
          m[i]=(float *) malloc(ncol*sizeof(float));
  
        /* return pointer to array of pointers to rows */
        return m;
}

void Profile(float *x, int ndat, float *y, float xcen, int xmin, 
		int xmax, float sigma)
{ 
	int i,k;
	float diff, denom, frac;

	denom = 1.0/sqrt(6.2832 * sigma * sigma);

/*		Below is denominator fro x^3	 */

/*	denom = 1.0/0.89298 * sigma * 2.0;	 */

	for (i=xmin,k=0; i<=xmax; i++, k++) {
	  y[k] = 0.0;
	  if(i >= 0 && i < ndat) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = fabs(xcen - x[i] + frac)/sigma;
/*
			Gaussian                     
	        y[k] += exp(-diff*diff/2.0 + 6.9077553)*denom; */
	        y[k] += exp(-diff*diff/2.0)*denom;
						

/*	        y[k] += exp(-diff*diff*diff)*denom;	*/

                }
	     y[k] /= 5.0;
	   }
	}
}

void Profilex2(float *x, int ndat, float *y, float xcen, int xmin, 
		int xmax, float sigma)
{ 
	int i,k;
	float diff, denom, frac;

	denom = 1.0/sqrt(6.2832 * sigma * sigma);

/*		Below is denominator for x^3	*/

/*	denom = 1.0/0.89298 * sigma * 2.0;	*/

	for (i=xmin,k=0; i<=xmax; i++, k++) {
	  y[k] = 0.0;
	  if(i >= 0 && i < ndat) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = fabs(xcen - x[i] + frac)/sigma;
/*
			Gaussian                     
	        y[k] += exp(-diff*diff/2.0 + 6.9077553)*denom; */
	        y[k] += diff*diff*exp(-diff*diff/2.0)*denom;
						

/*	        y[k] += exp(-diff*diff*diff)*denom;	*/

                }
	     y[k] /= 5.0;
	   }
	}
}

	
void FillXs(int *xmin, int *xmax,float *x,float *xcen, int ndat, int nTrace,
		float diff)
{
/*
//        Find max and min of profile influence
//	*/
   int current, place;

   current = 0;
   place = 0;

   while (current < nTrace) {
      while(x[place] < xcen[current] - diff && place < ndat) place++;
      if(place >= ndat) place--;
      xmin[current] = place;

      while(x[place] <= xcen[current] + diff && place < ndat) place++;
      if(place >= ndat) place--;
      xmax[current] = place - 1; 

      if(xmax[current] <= xmin[current]) {
         fprintf(stderr," Fiber # %d has no influence \n", current);
/*         ia[current] = 0;
         a[current] = 0.0; */
      } 

      place = xmin[current];
      current++;
   }     

}

void FillXs2(int *xmin, int *xmax,float *x,float *xcen, int ndat, int nTrace,
		float* sigma, float **abig)
{
/*
//        Find max and min of profile influence
//	*/
   int current, place, place2;
   float diff;

   current = 0;
   place = 0;

   while (current < nTrace) {
      diff = 4.0*sigma[current];
      while(x[place] < xcen[current] - diff && place < ndat) place++;
      if(place >= ndat) place--;

      place2 = place;
      while(x[place2] <= xcen[current] + diff && place2 < ndat) place2++;
      if(place2 >= ndat) place2--;

/*      if(place != xmin[current] || xmax[current] != place2-1) { 
	  if(abig[current*2] != 0x0) {
	       free(abig[current*2]);
	       abig[current*2] = 0x0;
	  }

	  if(abig[current*2+1] != 0x0) {
	       free(abig[current*2+1]);
	       abig[current*2+1] = 0x0;
	  }
	} */
      xmin[current] = place;
      xmax[current] = place2 - 1; 

      if(xmax[current] <= xmin[current]) {
         fprintf(stderr," Fiber # %d has no influence %f \n", current, diff);
/*         ia[current] = 0;
         a[current] = 0.0; */
      } 

      place = xmin[current];
      current++;
   }     

}

void FillABig(float **abig, float *x, float *xcen, int *xmin, int *xmax, 
                 int ndat, int nTrace, int nPoly, float sigma)
{
	int i, j, length;
        float *atemp;
	float xmid = 1024.0;
	float denom;

	denom = 1.0/sqrt(6.2832 * sigma * sigma);

        atemp = (float *)malloc(nPoly * sizeof(float));   

	for (j=0; j<nTrace; j++) {
	    length = xmax[j] - xmin[j] + 1;
	    if(length > 0) {
               abig[j] = (float *)malloc(length * sizeof(float));
	       if (sigma < 0) 
	         Profilex2(x, ndat, abig[j], xcen[j], xmin[j], xmax[j],-sigma); 
	       else
	         Profile(x, ndat, abig[j], xcen[j], xmin[j], xmax[j], sigma); 
	    }
	}
 
	for(j=0; j<nPoly; j++) 
           abig[j+nTrace] = (float *)malloc(ndat * sizeof(float));

        for (i=0;i<ndat;i++) {
	   chebyshevFunc((x[i]-xmid)/xmid, atemp, nPoly);
	   for(j=0; j<nPoly; j++) {
	     abig[j+nTrace][i] = atemp[j];
           }
        }

	free(atemp);
}

void FillABig2(float **abig, float *x, float *xcen, int *xmin, int *xmax, 
                 int ndat, int nTrace, int nPoly, float *sigma)
{
	int i, j, length;
	float xmid = 1024.0;
	float denom;
	float *atemp;


        atemp = (float *)malloc(nPoly * sizeof(float));   

	for (j=0; j<nTrace; j++) {
	    length = xmax[j] - xmin[j] + 1;
	    if(length > 0 && sigma[j] > 0.0) {
	       denom = 1.0/sqrt(6.2832 * sigma[j] * sigma[j]);

               abig[j*2] = (float *)realloc(abig[j*2],length * sizeof(float));
	       Profile(x, ndat, abig[j*2], xcen[j], xmin[j], xmax[j],sigma[j]); 

               abig[j*2+1] = (float *)realloc(abig[j*2+1],length*sizeof(float));
	       Profilex2(x, ndat,abig[j*2+1],xcen[j],xmin[j],xmax[j],sigma[j]); 
	    }
	}
 
	for(j=2*nTrace; j<2*nTrace + nPoly; j++) 
           abig[j] = (float *)realloc(abig[j], ndat * sizeof(float));

        for (i=0;i<ndat;i++) {
	   chebyshevFunc((x[i]-xmid)/xmid, atemp, nPoly);
	   for(j=0; j<nPoly; j++) {
	     abig[j+2*nTrace][i] = atemp[j];
           }
        }

	free(atemp);
}

/* Custom2 replaces lower triangle of a with sqrt(covar)  
	And we only loop over enough to work on 2 parameters for
	adjacent fibers   
   a is modified and cannot be used again for subsequent calls */
void cholslCustom2(float **a, int *ia, int nTrace, int nPoly, 
	                  float p[], float b[], float x[])
{
	int i,k;
	float sum;

	for(i=0;i<nTrace;i++) 
           if (ia[i]) {
	      sum = b[i];
	   
	      for (sum=b[i],k=i-1;k>=0 && k>=i-3;k--) 
                 if (ia[k])  sum -= a[i][k]*x[k];
	      x[i] = sum/p[i];
	}
	for(i=nTrace;i<nTrace+nPoly;i++) 
           if (ia[i]) {
	      for (sum=b[i],k=i-1;k>=0;k--) 
              if (ia[k])  sum -= a[i][k]*x[k];
	      x[i] = sum/p[i];
	}


	for(i=nTrace + nPoly - 1;i>=nTrace;i--) 
           if (ia[i]) {
	      for (sum=x[i],k=i+1;k<nTrace+nPoly;k++) 
                 if (ia[k])  sum -= a[k][i]*x[k];
	   x[i] = sum/p[i];
	}
	for(i=nTrace-1;i>=0;i--) 
           if (ia[i]) {
	      for (sum=x[i],k=i+1;k<nTrace && k<=i+3;k++) 
                 if (ia[k])  sum -= a[k][i]*x[k];
	      for (k=nTrace;k<nTrace + nPoly;k++) 
                 if (ia[k]) sum -= a[k][i]*x[k];
	   x[i] = sum/p[i];
	}

}

void cholslCustomCovar(float **a, int *ia, int nTrace, int nPoly, float p[])
{
	int i,k;
	float sum;
	float *x;
	int pl;
	
	for(pl=0;pl < nTrace+nPoly; pl++)  
           if (ia[pl]) {
	 
	      x = a[pl];

	      for(i=pl,sum=1.0;i<nTrace;i++,sum=0.0) 
                 if (ia[i]) {
	         for (k=i-1;k>=pl && k>=i-3;k--) 
                    if (ia[k]) sum-=a[i][k]*x[k];
	         x[i] = sum/p[i];
	      }
	   for(;i<nTrace+nPoly;i++,sum=0.0) 
              if (ia[i]) {
                 for (k=i-1;k>=pl;k--) 
                    if (ia[k]) sum -= a[i][k]*x[k];
	      x[i] = sum/p[i];
	   }

	   for(i=nTrace+nPoly-1; i>=nTrace && i >= pl;i--) 
              if (ia[i]) {
                 sum = x[i];
	         for (k=i+1;k<nTrace+nPoly;k++) 
                    if (ia[k]) sum -= a[k][i]*x[k];
	      x[i] = sum/p[i];
	   }
	   for(i=nTrace-1;i>=0 && i >= pl;i--) 
              if (ia[i]) {
                 sum = x[i];
	         for (k=i+1;k<nTrace && k<=i+3;k++) 
                    if (ia[k]) sum -= a[k][i]*x[k];
	         for (k=nTrace;k<nTrace + nPoly;k++) 
                    if (ia[k]) sum -= a[k][i]*x[k];
	         x[i] = sum/p[i];
	   }
	}

}

/*	Custom2 requires two (2) parameters per Trace in top left of covar */           
int choldcCustom2(float **a, int *ia, int nTrace, int nPoly, float p[])
{
	int i,j,k;
	float sum;

	for(i=0;i<nTrace;i++) 
	   if(ia[i]) {
	      for(j=i; j<nTrace && j<i+4; j++) 
	         if(ia[j]) {
                    for(sum = a[i][j], k=i-1; k >= 0 && k >= i-3; k--) 
	               if(ia[k]) sum -= a[i][k]*a[j][k];
              
	       if(i==j) {
	         if (sum <= 0.0) {
	           fprintf(stderr,"choldc failed %d %f\n",i,sum);
                   sum = a[i][j];
	           fprintf(stderr,"%f\n",sum);
                   for(k=i-1; k >= 0 && k >= i-3; k--) {
	              if(ia[k]) sum -= a[i][k]*a[j][k];
	              fprintf(stderr,"%f %f %f\n",sum, a[i][k], a[j][k]);
	           }
	          return -i;
	         }

	         p[i] = sqrt(sum);
	       } 
               else a[j][i] = sum/p[i];
	    }
	    for(j=nTrace; j<nTrace +nPoly; j++) 
	       if(ia[j]) {
                  for(sum = a[i][j], k=i-1; k >= 0 && k >= i-3; k--) 
	             if(ia[k]) sum -= a[i][k]*a[j][k];
              
                  a[j][i] = sum/p[i];
	    }

        }
	for(i=nTrace;i<nTrace+nPoly;i++) 
	   if(ia[i]) {
	   for(j=i;j<nTrace+nPoly;j++) 
	      if(ia[j]) {
	         for (sum=a[i][j],k=i-1;k>=0;k--) 
	            if(ia[k]) sum -= a[i][k]*a[j][k];
	         if(i==j) {
	            if (sum <= 0.0) fprintf(stderr,"choldc failed\n");
	            p[i] = sqrt(sum);
	      } else a[j][i] = sum/p[i];
           }
        }
	return 0;
}
           


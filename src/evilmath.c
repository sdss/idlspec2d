
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
   float    *  imerr,
   float       radius,
   float    *  xcen,
   float    *  xerr)
{
   long        ix;
   long        ix1;
   long        ix2;
   float       x1;
   float       x2;
   float       sumfx;
   float       sumf;
   float       sumss;
   float    *  window;

   /* Only recenter if the guess value is within the bounds of the image */
   if (*xcen > 0.0 && *xcen < nx) {

      /* Determine which pixel numbers over which to sum */
      x1 = *xcen - radius + 0.5;
      x2 = *xcen + radius + 0.5;
      ix1 = floor(x1);
      ix2 = floor(x2);

      /* If either end of the window is out of bounds in the image, then
         shrink both sides of the summing window until it is in bounds.
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
      nsum = ix2 - ix1 + 1;

      window = (float *) malloc(sizeof(float) *nsum);
      for (ix=ix1; ix <= ix2; ix++) {
         /* Determine the weight of a boxcar window function for this
          * pixel.  Note that the values of "window" will sum to 2*radius
          * unless the edge of the image is reached.
          */
         if (ix == ix1) {
            window[ix-ix1] = 1.0 + ix1 - x1;
         } else if (ix == ix2) {
            window[ix-ix1] = x2 - ix2;
         } else {
            window[ix-ix1] = 1.0;
         }
      }

      /* Compute the flux-weighted center. */
      sumfx = 0.0;
      sumf = 0.0;
      sumss = 0.0;
      for (ix=ix1; ix <= ix2; ix++) {
         sumfx += window[ix-ix1] * imrow[ix] * (ix - *xcen);
         sumf += window[ix-ix1] * imrow[ix];
         sumss += window[ix-ix1] * imerr[ix] * imerr[ix];
      }

      if (sumf > 0.0) {
         *xcen = sumfx / sumf + *xcen;
         *xerr = 0.0; /* NOT IMPLEMENTED ??? */
      } else {
         *xcen = -1.0;
         *xerr = 0.0;
      }

      free(window);
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

//		Below is denominator fro x^3

//	denom = 1.0/0.89298 * sigma * 2.0;

	for (i=xmin,k=0; i<=xmax; i++, k++) {
	  y[k] = 0.0;
	  if(i >= 0 && i < ndat) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = rabs(xcen - x[i] + frac)/sigma;
/*
			Gaussian                     
	        y[k] += exp(-diff*diff/2.0 + 6.9077553)*denom; */
	        y[k] += exp(-diff*diff/2.0)*denom;
						

//	        y[k] += exp(-diff*diff*diff)*denom;

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

//		Below is denominator fro x^3

//	denom = 1.0/0.89298 * sigma * 2.0;

	for (i=xmin,k=0; i<=xmax; i++, k++) {
	  y[k] = 0.0;
	  if(i >= 0 && i < ndat) {
	     for(frac = -0.4; frac <= 0.5; frac += 0.2)  {
	        diff = rabs(xcen - x[i] + frac)/sigma;
/*
			Gaussian                     
	        y[k] += exp(-diff*diff/2.0 + 6.9077553)*denom; */
	        y[k] += diff*diff*exp(-diff*diff/2.0)*denom;
						

//	        y[k] += exp(-diff*diff*diff)*denom;

                }
	     y[k] /= 5.0;
	   }
	}
}

	
void FillXs(int *xmin, int *xmax,float *x,float *xcen, int ndat, int nTrace,
		float diff)
{
//
//        Find max and min of profile influence
//
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
//
//        Find max and min of profile influence
//
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
           

int fixedGauss2(float x[], float y[], float ymod[], float invvar[],float *xcen, 
	int ndat, int nTrace, int nPoly, float *sigma, 
	float a[], int ia[], int ma, float **covar, float *chisq,
	int *xmin, int *xmax, float **abig)
{
/*     Stolen from NR lfit, but customized to fit fixed width gaussian's +
	a background chebyshev :
	Question, perform rejection inside here? */
/*	New version: 7/19/99, perform fit with fixed gaussian + a term
		for peturbed sigma */

        int i,j,k,l,m,mfit=0;
	int mfitTrace=0;
	int mfitPoly=0;
        float wt;
	int dTrace = 2*nTrace;
	int lm, extra, extra2;
        float *ysub = (float *)malloc(ndat * sizeof(float));   
        float *p= (float *)malloc(ma * sizeof(float));   
        float *beta = (float *)malloc(ma * sizeof(float));   
        float *a2= (float *)malloc(ma * sizeof(float));   

//
//	  Find basis functions (stored in ABig)
//	  Find max and min of profile influence
//	  Fill covar and beta matrices
//

/*  Expect xcen and x to be monotonically increasing  */


      FillXs2(xmin, xmax, x, xcen, ndat, nTrace, sigma, abig);

	printf("hello\n");
      FillABig2(abig, x, xcen, xmin, xmax, ndat, nTrace, nPoly, sigma);
	printf("hello2\n");

      CheckFibers(abig, xmin, xmax, nTrace, a, ia, invvar);
	printf("hello3\n");

        for (j=0;j<ma;j++) {
                if (ia[j]) mfit++;
	        }
        if (mfit == 0) {
	   printf("lfit: no parameters to be fitted\n");
	   return 3;
	}

        printf("mfit: %d\n",mfit);

        for (j=0;j<dTrace;j++) if (ia[j]) mfitTrace++;
        for (j=dTrace;j<dTrace+nPoly;j++) if (ia[j]) mfitPoly++;

//  zero out entire covar matrix

        for (j=0;j<ma;j++) {
                for (k=0;k<ma;k++) covar[j][k]=0.0;
                beta[j]=0.0;
        }

/* Subtract out fixed variables first */
     for(i=0;i<ndat;i++) {
         ysub[i] = y[i];
     }

     if(mfit < ma) {	
        for (j=0;j<dTrace;j++) 
           if (!ia[j]) 
              for (i=xmin[j/2],k=0;i<=xmax[j/2];i++,k++) 
                 ysub[i] -= a[j] * abig[j][k];
        
        for (j=dTrace;j<dTrace+nPoly;j++) 
           if (!ia[j]) 
              for (i=0;i<ndat;i++) 
                 ysub[i] -= a[j] * abig[j][i];
     }
           
   
/* Fill dTrace x dTrace  of covar first  */

  for (l=0,j=0;l<nTrace;l++) 
    for (extra=0; extra < 2; extra++)
       j=2*l+extra;
       if (ia[j]) {
         for (i=xmin[l];i<=xmax[l];i++) {
            wt = abig[j][i-xmin[l]] * invvar[i];
            beta[j] += wt * ysub[i];
            covar[j][j] += wt * abig[j][i-xmin[l]];

/*	Overlap with with chebyshev polynomial*/

            for (m = dTrace+nPoly-1; m>=dTrace; m--) 
                if(ia[m]) covar[j][m] += wt*abig[m][i];
         }

/*	Fill in symmetric pieces   */

         for (m = dTrace+nPoly-1; m>=dTrace; m--) 
                if(ia[m]) covar[m][j] = covar[j][m];
	

/*	Overlap with previous fiber  */

	 if (j > 0) {
	    k= j;
	    for (extra2 = extra; extra2 > -2; extra2--) 
              k=2*l+extra2-1;
              if (ia[k])  {
	       lm = l - 1;
	       if(extra2 > 0) lm = l;
	       if(lm >= 0) {
                  for (i = xmin[l];i<=xmax[lm];i++) {
                     wt = abig[j][i-xmin[l]] * invvar[i];
                     covar[j][k] += wt * abig[k][i-xmin[lm]];
	          }
                  covar[k][j] = covar[j][k];
               }
	    }
	 }

      }
/*   printf("nTrace done\n");

 Fill nPoly x nPoly of covar */
   for (l = dTrace+nPoly-1 , j=mfit; l >= dTrace ; l--) 
      if (ia[l]) {
         j--;
         for (i=0; i<ndat; i++) {
            wt = abig[l][i] * invvar[i];
            beta[j] += wt * ysub[i];
            for (k=mfit-1,m=dTrace+nPoly-1; m>=l; m--)
               if (ia[m]) covar[j][k--] += wt * abig[m][i];
         }
         for (k=mfit-1,m=dTrace+nPoly-1; m>l; m--)
            if (ia[m]) covar[k][j] = covar[j][k--];
      }
   printf("nPoly done\n");
   
   choldcCustom2(covar, ia, dTrace, nPoly, p); 
   printf("choldc Custom2 done\n");

   cholslCustom2(covar, ia, dTrace, nPoly, p, beta, a2); 
   printf("cholsl done\n");

   cholslCustomCovar(covar, ia, dTrace, nPoly, p); 
   printf("cholsl Covar done\n");
         
//   printf("Gaussj done\n");

   for (j=-1,l=0;l<ma;l++)
       if (ia[l]) {
	   a[l]=a2[++j];
//	   a[l]=beta[++j];
//           printf("%d %f %f \n", l, a[l], a2[l]);
//	   error[l] = 1.0/p[j];
       }


   *chisq=0.0;

   for(i=0;i<ndat;i++) ymod[i] = 0.0;

   for (l=0;l<dTrace;l++) {
//      fprintf(stderr,"%d %f %f\n", l, a[l], error[l]);
      for (i=xmin[l/2], k=0; i <= xmax[l/2]; i++, k++) 
         ymod[i] += a[l]*abig[l][k];
    }

   for (l=dTrace;l<dTrace+nPoly;l++) {
//      fprintf(stderr,"%d %f %f\n", l, a[l], error[l]);
      for(i=0; i < ndat; i++) 
         ymod[i] += a[l]*abig[l][i];
    }

   for (i=0;i<ndat;i++) *chisq += (y[i]-ymod[i])*(y[i]-ymod[i])*invvar[i];

   printf("Chisq: %f\n", *chisq);

	free(beta);
	free(p);
	free(ysub);
	free(a2);
	printf("wrapping up\n");
        return 1;
}

void chebyshevFunc(float x, float *coeff, int nCoeff)
{
	/* Return Terms of Chebyshev Polynomial */

        int i;
        float twox;

        if(nCoeff > 0) coeff[0] = 1.0;
        if(nCoeff > 1) coeff[1] = x;

        twox = 2.0*x;
        for(i=2;i<nCoeff;i++)
           coeff[i] = twox * coeff[i-1] - coeff[i-2];

      return;
}



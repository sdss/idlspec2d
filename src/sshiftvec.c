#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"

IDL_LONG sshiftvec
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float    *  shiftvec;
   float    ** simage;
   IDL_LONG    sincrad;
   float       dampfac;
   float       eps;

   long        ix;
   long        iy;
   long        ip;
   long        ishift;
   long        xstart;
   long        xend;
   long        sincsize;
   long        nleft;
   long        ix1;
   long        ix2;
   long        nright;
   float       fshift;
   float       y;
   float       py;
   float       dampfac2;
   float    *  sincvec;
   IDL_LONG    retval = 1;
   const float PI = 3.141592654;

   /* Allocate pointers from IDL */
   nx = *((IDL_LONG *)argv[0]);
   ny = *((IDL_LONG *)argv[1]);
   fimage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) fimage[iy] = (float *)argv[2] + iy*nx;
   shiftvec = (float *)argv[3];
   simage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) simage[iy] = (float *)argv[4] + iy*nx;
   sincrad = *((IDL_LONG *)argv[5]);
   dampfac = *((float *)argv[6]);
   eps = *((float *)argv[7]);

   dampfac2 = dampfac * dampfac;

   /* Allocate memory for the sinc vector */
   sincsize = 2 * sincrad + 1;
   sincvec = malloc(sincsize * sizeof(float));

   /* Loop through each row of the image */
   for (iy=0; iy < ny; iy++) {

      /* First, split the desired shift into a fractional and integer part. */
      ishift = shiftvec[iy]; /* Convert from float to int */
      fshift = shiftvec[iy] - ishift;

      /* Determine the valid limits in the output array */
      if (ishift > 0) {
         xstart = ishift;
         xend = nx-1;
      } else {
         xstart = 0;
         xend = nx-1+ishift;
      }

      if ( fabs(fshift) > eps && fabs(fshift) < 1.0-eps ) { /* FRACTIONAL SHIFT */

         /* Initialize the sinc vector */
         for (ip=0; ip < sincsize; ip++) {
            y = fshift - (ip - sincrad);
            py = PI * y;
            sincvec[sincsize-1-ip] = exp(-y*y/dampfac2) * sin(py)/py;
         }

         /* Loop through each column of the image and convolve the sinc array
          * with the input data.  This is the shift. */
         for (ix=xstart; ix <= xend; ix++) {

            simage[iy][ix] = 0.0;

            /* [ix1,ix2] describe the X limits in the input array */
            ix1 = ix - ishift - sincrad;
            ix2 = ix1 + sincsize;

            nleft = (ix1 < 0) ? -ix1 : 0;
            nright = (ix2 >= nx) ? ix2 - nx + 1 : 0;
            for (ip=0; ip < nleft; ip++) {
               simage[iy][ix] += fimage[iy][0] * sincvec[ip];
            }
            for (ip=nleft; ip < sincsize-nright; ip++) {
               simage[iy][ix] += fimage[iy][ix1+ip] * sincvec[ip];
            }
            for (ip=sincsize-nright; ip < sincsize; ip++) {
               simage[iy][ix] += fimage[iy][nx-1] * sincvec[ip];
            }
         }

      } else { /* INTEGER SHIFT */

         for (ix=xstart; ix <= xend; ix++) {
            simage[iy][ix] = fimage[iy][ix-ishift];
         }

      }

      /* Copy values on the ends for out-of-bounds pixels */
      for (ix=0; ix < xstart; ix++) {
         simage[iy][ix] = simage[iy][xstart];
      }
      for (ix=xend+1; ix < nx; ix++) {
         simage[iy][ix] = simage[iy][xend];
      }
   }

   /* Free temporary memory */
   free(sincvec);
   free(fimage);
   free(simage);

   return retval;
}


#include <stdio.h>
#include <math.h>
#include "export.h"
#include "evilmath.h"

/* Recenter one position in a row of data using a flux-weighted center */
void recenter_fweight
  (IDL_LONG    nx,
   float    *  imrow,
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
   float       window;

   /* Only recenter if the guess value is within the bounds of the image */
   if (*xcen > 0.0 && *xcen < nx) {

      /* Determine which column numbers over which to sum */
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

      /* Compute the flux-weighted center. */
      sumfx = 0.0;
      sumf = 0.0;
      for (ix=ix1; ix <= ix2; ix++) {
         /* Determine the weight of a boxcar window function for this
          * pixel.  Note that the values of "window" will sum to 2*radius
          * unless the edge of the image is reached.
          */
         if (ix == ix1) {
            window = 1.0 + ix1 - x1;
         } else if (ix == ix2) {
            window = x2 - ix2;
         } else {
            window = 1.0;
         }
         sumfx += window * imrow[ix] * (ix - *xcen);
         sumf += window * imrow[ix];
      }

      if (sumf > 0.0) {
         *xcen = sumfx / sumf + *xcen;
         *xerr = 0.0; /* NOT IMPLEMENTED ??? */
      } else {
         *xcen = -1.0;
         *xerr = 0.0;
      }
   }
}


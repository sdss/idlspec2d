#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"

IDL_LONG extract_boxcar
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float       radius;
   IDL_LONG    ncen;
   float    *  xcen;
   IDL_LONG *  ycen;
   float    *  fextract;

   long        ix;
   long        ix1;
   long        ix2;
   long        iy;
   long        icen;
   float       x1;
   float       x2;
   float       window;
   IDL_LONG    retval = 1;

   /* Allocate pointers from IDL */
   nx = *((IDL_LONG *)argv[0]);
   ny = *((IDL_LONG *)argv[1]);
   fimage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) fimage[iy] = (float *)argv[2] + iy*nx;
   radius = *((float *)argv[3]);
   ncen = *((IDL_LONG *)argv[4]);
   xcen = ((float *)argv[5]);
   ycen = ((IDL_LONG *)argv[6]);
   fextract = ((float *)argv[7]);

   /* Loop through each center value */
   for (icen=0; icen < ncen; icen ++) {

      /* Determine which column numbers over which to sum */
      x1 = xcen[icen] - radius + 0.5;
      x2 = xcen[icen] + radius + 0.5;
      ix1 = floor(x1);
      ix2 = floor(x2);

      /* Only extract if some of the pixels are within bounds of the image */
      if (x1 < nx-1 && x2 > 0.0) {

         /* If either end of the window is out of bounds in the image, then
          * truncate that side of the extraction window at the image boundary.
          */
         if (x1 < 0.0) {
            x1 = 0.0;
            ix1 = 0;
         }
         if (x2 > nx-1) {
            x2 = nx-1;
            ix2 = nx-1;
         }

         /* Extract */
         fextract[icen] = 0.0;
         for (ix=ix1; ix <= ix2; ix++) {
            /* Determine the weight of a boxcar window function for
             * this pixel.
             */
            if (ix == ix1 && ix == ix2) {
               window = x2 - x1;
            } else if (ix == ix1) {
               window = 1.0 + ix1 - x1;
            } else if (ix == ix2) {
               window = x2 - ix2;
            } else {
               window = 1.0;
            }
            fextract[icen] += window * fimage[ycen[icen]][ix];
         }
      }
   }

   /* Free temporary memory */
   free(fimage);

   return retval;
}

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "evilmath.h"

IDL_LONG trace_crude
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   IDL_LONG    ny;
   float    ** fimage;
   float       radius;
   IDL_LONG    ntrace;
   float    *  xstart;
   IDL_LONG *  ystart;
   float    ** xset;

   long        iy;
   long        itrace;
   float       xerr;
   IDL_LONG    retval = 1;

   /* Allocate pointers from IDL */
   nx = *((IDL_LONG *)argv[0]);
   ny = *((IDL_LONG *)argv[1]);
   fimage = (float **)malloc(ny * sizeof(float *)); /* build pointers only */
   for (iy=0; iy < ny; iy++) fimage[iy] = (float *)argv[2] + iy*nx;
   radius = *((float *)argv[3]);
   ntrace = *((IDL_LONG *)argv[4]);
   xstart = ((float *)argv[5]);
   ystart = ((IDL_LONG *)argv[6]);
   xset = (float **)malloc(ntrace * sizeof(float *)); /* build pointers only */
   for (itrace=0; itrace < ntrace; itrace++)
    xset[itrace] = (float *)argv[7] + itrace*ny;

   /* Loop through each trace */
   for (itrace=0; itrace < ntrace; itrace ++) {

      /* RECENTER INITIAL ROW */
      iy = ystart[itrace];
      xset[itrace][iy] = xstart[itrace];
      recenter_fweight(nx, fimage[iy], radius, &xset[itrace][iy], &xerr);

      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO LARGER ROW NUMBERS */
      for (iy=ystart[itrace]+1; iy < ny; iy++) {
         xset[itrace][iy] = xset[itrace][iy-1];
         recenter_fweight(nx, fimage[iy], radius, &xset[itrace][iy], &xerr);
      }

      /* LOOP FROM INITIAL (COL,ROW) NUMBER TO SMALLER ROW NUMBERS */
      for (iy=ystart[itrace]-1; iy >= 0; iy--) {
         xset[itrace][iy] = xset[itrace][iy+1];
         recenter_fweight(nx, fimage[iy], radius, &xset[itrace][iy], &xerr);
      }

   }

   /* Free temporary memory */
   free(fimage);
   free(xset);

   return retval;
}


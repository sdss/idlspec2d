#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "evilmath.h"
  

IDL_LONG extract_row
  (int      argc,
   void *   argv[])
{
   IDL_LONG    nx;
   float     * xdummy;
   float     * fimage;
   float     * invvar;
   float     * ymod;
   IDL_LONG    nTrace;
   IDL_LONG    nPoly;
   IDL_LONG    nCoeff;
   IDL_LONG    ma;
   float     * xcen;
   float     * sigma;
   IDL_LONG  * ia;
   float     * ans;
   float     * p;
   float    ** covar;
   float     * fscat;

   long        iy;
   IDL_LONG    retval = 1;

   int	       argct;

   /* Allocate pointers from IDL */

   argct  = 0;
   nx     = *((IDL_LONG *)argv[argct++]);
   xdummy = (float *)argv[argct++];
   fimage = (float *)argv[argct++];
   invvar = (float *)argv[argct++];
   ymod   = (float *)argv[argct++];

   nTrace = *((IDL_LONG *)argv[argct++]);
   nPoly  = *((IDL_LONG *)argv[argct++]);

   xcen   = (float *)argv[argct++];
   sigma  = (float *)argv[argct++];

   nCoeff = *((IDL_LONG *)argv[argct++]);
   ma     = *((IDL_LONG *)argv[argct++]);
   ans    = (float *)argv[argct++];
   ia     = (IDL_LONG *)argv[argct++];
   p      = (float *)argv[argct++];
   fscat  = (float *)argv[argct++];

   covar  = (float **)malloc(ma * sizeof(float *)); 
   for (iy=0; iy < ma; iy++) covar[iy] = (float *)argv[argct]+iy*ma;

   fprintf(stderr, "Going to fit_row\n");

   /* Free temporary memory */
   free(covar);

   return retval;
}


;------------------------------------------------------------------------------
;+
; NAME:
;   locateskylines_test
;
; PURPOSE:
;   Compute shifts for arc lines used in wavelength calibration.
;
; CALLING SEQUENCE:
;   locateskylines, skylinefile, fimage, ivar, wset, xarc, arcshift=arcshift,
;    [ nord=, xsky=, ysky=, skywaves=, skyshift= ]
;
; INPUTS:
;   skylinefile - Name of skyline file (with air wavelengths in Angstroms)
;   fimage      - Flattened image containing skylines [NPIX,NFIBER]
;   ivar        - Inverse variance of FIMAGE
;   wset        - Wavelength solution traceset (pix -> lambda)
;   xarc        - Arc line positions ???
;
; OPTIONAL INPUTS:
;   nord     - Order of fit to delta x positions; default to 4
;
; OUTPUTS:
;   arcshift    - Shifts to apply to arc lines in pix [NROW,NTRACE]
;
; OPTIONAL OUTPUTS:
;   xsky        - Pixel position of sky lines [NFIBER,NLINE]
;   skywaves    - Wavelengths of sky lines used for computing shifts (Ang)
;   skyshift    - Shifts to apply to sky lines in pix [NROW,NTRACE]
;
; COMMENTS:
;   The wavelength as a function of fiber number is only allowed to
;   vary quadratically.  The scale is only allowed to vary by the same
;   factor for the entire image.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   trace_fweight()
;   trace_gweight()
;   traceset2pix()
;   traceset2xy
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO
;   18-Nov-1999  Moved skyline QA to fit_skyset (SMB)
;-
;------------------------------------------------------------------------------

pro locateskylines_test, skylinefile, fimage, ivar, wset, xarc, arcshift=arcshift, $
 xsky=xsky, skywaves=skywaves, skyshift=skyshift

   if (NOT keyword_set(nord)) then nord = 4

   splog, 'Reading sky line file ', skylinefile
   readcol, skylinefile, skywaves, /silent, form='d'

   npix = (size(fimage))[1]
   nfiber = (size(fimage))[2]
   ncoeff = (size(wset.coeff))[1]

   ;---------------------------------------------------------------------------
   ; Make xsky, ysky, pairs from skywaves
   ;---------------------------------------------------------------------------
   ; xpredict contains skyline positions based on arc fit

   nskyline = n_elements(skywaves)
   ysky = dindgen(nfiber) # (dblarr(nskyline)+1)
   xpredict = transpose( traceset2pix(wset, alog10(skywaves)) )

   ;---------------------------------------------------------------------------
   ; Discard lines that are out of bounds for any fiber
   ;---------------------------------------------------------------------------

   qgood = total(((xpredict LE 0) OR (xpredict GE npix)),1) EQ 0 
   gind = where(qgood, nskyline)

   if (nskyline EQ 0) then begin
      splog, 'WARNING: No sky lines in wavelength range'
      return
   endif

   ; Trim list
   xpredict = xpredict[*,gind]
   ysky = ysky[*,gind]
   skywaves = skywaves[gind]

   ;---------------------------------------------------------------------------
   ; Trace the sky lines
   ;---------------------------------------------------------------------------

   ;----------
   ; Recenter on every sky line in every fiber, using the predicted positions
   ; from the wavelength solution.

   xskytmp = trace_fweight(fimage, xpredict, ysky, invvar=ivar, $
    xerr=xerr, radius=3.0)

   ;----------
   ; Apply a global shift, and iterate this re-centering twice.

   for iiter=0, 1 do begin
      medianshift = median(xskytmp-xpredict)
      xskytmp = trace_fweight(fimage, xpredict+medianshift, $
                              ysky, invvar=ivar, radius=3.0, xerr=fxerr)
   endfor

   ;----------
   ; Recenter using a gaussian fit.

   xsky = trace_gweight(fimage, xskytmp, ysky, invvar=ivar, sigma=1.0, $
    xerr=gxerr) 

   ;---------------------------------------------------------------------------
   ; Discard bad pixels
   ;---------------------------------------------------------------------------

   ;----------
   ; Create mask of centroids with large errors.
   ; Reject if errors from flux-weight-centroids or gaussian centroids
   ; are large, or if these two methods disagree.

maxerr = 0.05 ; ???
   mask = gxerr LT maxerr AND gxerr LT maxerr AND abs(xsky-xskytmp) LT maxerr
   igood = where(mask, ngood)
   if (ngood LT 100) then begin
      splog, 'WARNING: Too few good sky centroids for shifting wavelengths'
      return
   endif

   ;----------
   ; Determine how many of the sky lines have good S/N
   ; in at least 50% of the fibers

   junk = where( total(mask,1) GT 0.5 * nfiber, ngoodsn)

   if (ngoodsn EQ 0) then begin
      splog, 'WARNING: No sky lines with good S/N'
      return
   endif

   ;---------------------------------------------------------------------------
   ; Fit the shifts
   ;---------------------------------------------------------------------------

   ;----------
   ; First fit a quadratic shift, as a function of fiber number, to the
   ; offsets

   fibernums = djs_laxisgen([nfiber,nskyline], iaxis=0)
   res1 = polyfitw(float(fibernums[igood]), (xsky-xpredict)[igood], $
    1./(gxerr[igood])^2, nord, xfit1)
;splot, xpredict[igood], (xsky-xpredict)[igood],ps=3
;soplot, xpredict[igood], xfit1, ps=3, color='red'

   ;----------
   ; Now fit for a scale change as a the d(scale)/d(x) on the CCD.
   ; Only fit this if there are at least 3 sky lines with good S/N.

   if (ngoodsn GE 3) then begin
      res2 = polyfitw(float(xpredict[igood]), (xsky-xpredict)[igood]-xfit1, $
       1./(gxerr[igood])^2, 1, xfit2)
;soplot, xpredict[igood], (xsky-xpredict)[igood]-xfit1,ps=3,color='red'
;soplot, xpredict[igood], xfit2
   endif else begin
      res2 = 0
   endelse

   ;---------------------------------------------------------------------------
   ; Compute shift for sky line positions.

   skyshift = poly(fibernums, res1) + poly(xpredict, res2) ; The sign ???

   ;---------------------------------------------------------------------------
   ; Compute shift for arc line positions.

   arcfibnums = djs_laxisgen(size(xarc,/dimens), iaxis=1)
   arcshift = poly(arcfibnums, res1) + poly(xarc, res2) ; The sign???

   splog, 'Median arc shift = ', median(arcshift), ' pix'
   djs_iterstat, arcshift, sigma=sig
   splog, 'Dispersion in arc shift = ', sig, ' pix'

   return
end

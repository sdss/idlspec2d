;------------------------------------------------------------------------------
;+
; NAME:
;   locateskylines
;
; PURPOSE:
;   Compute shifts for arc lines used in wavelength calibration.
;
; CALLING SEQUENCE:
;   locateskylines, skylinefile, fimage, ivar, wset, xarc, arcshift=arcshift,
;    [ xsky=, ysky=, skywaves=, skyshift= ]
;
; INPUTS:
;   skylinefile - Name of skyline file (with air wavelengths in Angstroms)
;   fimage      - Flattened image containing skylines [NPIX,NFIBER]
;   ivar        - Inverse variance of FIMAGE
;   wset        - Wavelength solution traceset (pix -> lambda)
;   xarc        - Arc line positions from arc frame [pix]
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
;   vary quadratically.  The scale 
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   fitmeanx()
;   trace_fweight()
;   trace_gweight()
;   traceset2pix()
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO
;   18-Nov-1999  Moved skyline QA to fit_skyset (SMB)
;-
;------------------------------------------------------------------------------

pro locateskylines, skylinefile, fimage, ivar, wset, xarc, arcshift=arcshift, $
 xsky=xsky, skywaves=skywaves, skyshift=skyshift

   splog, 'Reading sky line file ', skylinefile
   skywaves = 0
   readcol, skylinefile, skywaves, /silent, form='D', numline=numline
   if (NOT keyword_set(skywaves)) then $
    message, 'Unable to read sky wavelengths from file ' + skylinefile

   dims = size(fimage, /dimens)
   npix = dims[0]
   nfiber = dims[1]
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
      skyshift = 0
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
   ; Fit the shifts
   ;---------------------------------------------------------------------------

   lambda = alog10(skywaves)
   xskyold = xsky
   inmask = gxerr LT 900.
   ; Insist that 80% of the sky lines be well-fit to use that line at all
   for i=0, nskyline-1 do $
    inmask[*,i] *= mean(inmask[*,i]) GT 0.80
   xdiff = fitmeanx(wset, lambda, xskyold, aveinvvar, inmask=inmask, mx=mx)

   junk = where(aveinvvar GT 1./(0.2)^2, ngood)
   if (ngood EQ 0) then begin
      splog, 'WARNING: No good sky lines (with std-dev less than 0.2 pixels)'
      skyshift = 0
      return
   endif else if (ngood LE 3) then begin
      shiftcoeff = 1 ; Fit an offset for each fiber
   endif else begin
      shiftcoeff = 2 ; Fit an offset + scale factor for each fiber
   endelse

   splog, 'Fitting ', shiftcoeff, ' coefficients for wavelength shifts'

   invvar = 0.0 * xdiff
   for i=0, nskyline-1 do invvar[*,i] = aveinvvar[i]
   xy2traceset, transpose(mx), transpose(xdiff), shiftset, ncoeff=shiftcoeff, $
    invvar=transpose(invvar), xmin=0, xmax=npix-1, maxiter=0

   ;---------------------------------------------------------------------------
   ; Compute shift for sky line positions.

   traceset2xy, shiftset, transpose(xsky), skyshift
   skyshift = transpose(skyshift)

   ;---------------------------------------------------------------------------
   ; Compute shift for arc line positions.

   traceset2xy, shiftset, transpose(xarc), arcshift
   arcshift = transpose(arcshift)

   splog, 'Median arc shift = ', median(arcshift), ' pix'
   djs_iterstat, arcshift, sigma=sig
   splog, 'Dispersion in arc shift = ', sig, ' pix'

   maxshift = max(abs(arcshift))
   if (maxshift GT 3.5) then begin
      splog, 'WARNING: Maximum sky-line shift is ', maxshift, ' (DISABLING)'
      skyshift = 0
      arcshift = 0
   endif else begin
      splog, 'Maximum sky-line shift is ', maxshift
   endelse

   return
end

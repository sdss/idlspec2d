;------------------------------------------------------------------------------
;+
; NAME:
;   locateskylines
;
; PURPOSE:
;   Tweak wavelength calibration using skylines
;
; CALLING SEQUENCE:
;   locateskylines, skylinefile, fimage, ivar, wset, xsky, ysky, skywaves, $
;    lambda=lambda, xcoeff=xcoeff
;
; INPUTS:
;   skylinefile - filename of skyline file
;   fimage      - flattened image containing skylines (npix x nfiber)
;   ivar        - inverse variance of fimage
;   wset        - traceset (pix -> lambda)
;
; OPTIONAL KEYWORDS:
;   lambda      - override skyline file (log10 Angstroms)
;
; OUTPUTS:
;   xsky        - pixel position of lines [nfiber, nline]
;   ysky        - fiber number [nfiber, nline]
;   skywaves    - wavelengths used (Angstroms)
;
; OPTIONAL OUTPUTS:
;   errcode     - returns errcode (see below)
;
; COMMENTS:
;   Error codes:
;  -1 - <unknown error>
;   0 - Everything fine
;
; EXAMPLES:
;
; BUGS:
;  I've now changed traceset2pix to use the forward trace set, but have
;  not yet made any changes to this code.
;
; PROCEDURES CALLED:
;  fitmeanx         - internal routine for FITARCIMAGE
;  trace_fweight()
;  traceset2pix()
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO
;   18-Nov-1999  Moved skyline QA to fit_skyset (SMB)
;-
;------------------------------------------------------------------------------

pro locateskylines, skylinefile, fimage, ivar, wset, $
 xsky, ysky, skywaves, lambda=lambda, xcoeff=xcoeff
  

   if (keyword_set(lambda)) then begin 
      skywaves = 10.d^lambda
   endif else begin 
      print, 'Reading ', skylinefile
      readcol, skylinefile, skywaves, /silent, form='d'
   endelse

   nskyline = n_elements(skywaves)
   npix = (size(fimage))[1]
   nfiber = (size(fimage))[2]
   ncoeff = (size(wset.coeff))[1]

   ;---------------------------------------------------------------------------
   ; Make xsky, ysky, pairs from skywaves
   ;---------------------------------------------------------------------------
   ; xarc contains skyline positions based on arc fit

   ysky = dindgen(nfiber) # (dblarr(nskyline)+1)
   xarc = transpose( traceset2pix(wset, alog10(skywaves)) )

   ;---------------------------------------------------------------------------
   ; Discard lines that are out of bounds
   ;---------------------------------------------------------------------------
   good = total(((xarc LE 0) OR (xarc GE npix)),1) EQ 0 
   gind = where(good)
   if (gind[0] EQ -1) then message, 'No sky lines!!!'

   ; trim list
   xarc = xarc[*,gind]
   ysky = ysky[*,gind]
   skywaves = skywaves[gind]

   ; Iterate trace_fweight ???
   xskytmp = trace_fweight(fimage, xarc, ysky, invvar=ivar) 
   xsky = trace_fweight(fimage, xskytmp, ysky, invvar=ivar, radius=2.0) 

   lambda = alog10(skywaves)
   xskyold = xsky
   xdiff = fitmeanx(wset, alog10(skywaves), xskyold, aveinvvar, mx = mx )

   good = where(aveinvvar[0,*] GT 1.0, ngood)
   if (good[0] EQ -1) then begin
     splog, 'No good sky lines (with residuals less than 0.2 pixels'
     return
   endif

   shiftcoeff = 1
   if (ngood GT 10) then shiftcoeff = 2
   xcoeff = dblarr(nfiber,shiftcoeff+1)
   for i=0,nfiber - 1 do  $
     xcoeff[i,*] = polyfitw(mx[i,good],xdiff[i,good], aveinvvar[i,good], $
                 shiftcoeff)
   
   xsky = xsky[*,good] 
   skywaves = skywaves[good]

      
   return

end

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
;    lambda=lambda, errcode=errcode
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
;-
;------------------------------------------------------------------------------

pro locateskylines, skylinefile, fimage, ivar, wset, $
 xsky, ysky, skywaves, lambda=lambda, errcode=errcode

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

   fmed=fimage-fimage
   for i=0, nfiber-1 do fmed[*,i] = median(fimage[*,i], 17)

   ; Iterate trace_fweight ???
   xskytmp = trace_fweight(fimage, xarc, ysky, invvar=ivar) 
   xsky = trace_fweight(fimage, xskytmp, ysky, invvar=ivar) 

   ;---------------------------------------------------------------------------
   ;  Bonehead Statistics
   ;---------------------------------------------------------------------------

   ; xsky now contains positions of skylines.  We have no idea if they 
   ; are any good unless we do some checks.

   nskyline = (size(xsky))[2]
   mean  = fltarr(nskyline)
   sigma = fltarr(nskyline)
   xres  = xarc - xsky

   ; Amp 1
   for i=0, nskyline-1 do begin 
      djs_iterstat, xres[3:150,i], mean=mn, sigma=sig ; HORRIBLE HARDWIRE???
      mean[i] = mn
      sigma[i] = sig
   endfor 
   mean1 = mean
   sigma1 = sigma

   ; Amp 2
   for i=0, nskyline-1 do begin 
      djs_iterstat, xres[170:316,i], mean=mn, sigma=sig ; HORRIBLE HARDWIRE???
      mean[i] = mn
      sigma[i] = sig
   endfor 
   mean2 = mean
   sigma2 = sigma

   ; Both amps
   for i=0, nskyline-1 do begin 
      djs_iterstat, xres[*,i], mean=mn, sigma=sig
      mean[i] = mn
      sigma[i] = sig
   endfor 

   ;---------------------------------------------------------------------------
   ;  Quality assurance
   ;---------------------------------------------------------------------------

   pmulti = !p.multi
   !p.multi = [0,1,2]
   plot, skywaves, mean1, yr=[-.3,.3]+median(mean), /yst, $
    xtit='Wavelength (A)', ytit='Delta x (Pix)', /xst, $
    title='Sky line residual'
   errplot, skywaves, mean1-sigma1, mean1+sigma1
   oplot, skywaves, mean2, line=1

   ; Reject outliers with "drop dead" requirements
   good = (abs(mean) LT 3) AND (abs(sigma) LT 1.0)
   gind = where(good, nline)   
   if (gind[0] EQ -1) then message, 'No sky lines!!!'
   print,nline, ' good lines - ',fix(total(good EQ 0)), ' rejected'

   ; Trim list again 
   mean  = mean[gind]
   sigma = sigma[gind]
   xarc  = xarc[*,gind]
   xsky  = xsky[*,gind]
   ysky  = ysky[*,gind]
   skywaves = skywaves[gind]
      
   if (nline GE 1) then nord=0
   if (nline GE 3) then nord=1
   if (nline GE 7) then nord=2

   print, 'Fit order:', nord

   flexcoeff = polyfitw(skywaves, mean, 1./sigma^2, nord, yfit)
   oplot, skywaves, yfit, line=2

   infostr = string('Dispersion:', stddev(mean-yfit)*1E3,' mpix', $
    format='(A,F7.1,A)')

   plot, skywaves, mean-yfit, $
    xrange=[min(skywaves)-100,max(skywaves)+100], yr=[-.2,.2], $
    xtit='Wavelength (A)', ytit='Delta x (Pix)', $
    title='After flexure correction'
   errplot,skywaves,mean-yfit-sigma,mean-yfit+sigma
   xyouts, 0.95, 0., systime(), /normal, align=1, chars=0.5
   xyouts, 0.05, 0., infostr, /norm
   print, infostr

print, 'No flexure correction...yet'

;
;   SMB (10/31/99): Prepare for flexure correction here
;   Return the best skyline positions, and fit in spreduce
;
   lambda = alog10(skywaves)
   xskyold = xsky
   xsky = fitmeanx(wset, lambda, xskyold)

   !p.multi=pmulti
   return

end

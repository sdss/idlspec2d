function quickwave, arcname, flatname, outarc, radius=radius

   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(outarc) NE 1) then return, 0
   if (n_elements(flatname) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0

   ;----------
   ; Read in image

   sdssproc, arcname, arcimg, hdr=hdr, color=color
   camname = strtrim(sxpar(hdr, 'CAMERAS'),2)

   ;----------
   ; Read in the reduced flat

   tset = mrdfits(flatname,2)
   fibermask = mrdfits(flatname,4)
   traceset2xy, tset, ycen, xcen 

   ;----------
   ; Boxcar extract the arc

   flux = extract_boxcar(arcimg, xcen, radius=radius)

   ;----------
   ; Estimate inverse variance

   fluxivar = 1.0 / (abs(flux) + 10.0)

   ;----------
   ; Now find the wavelength solution

   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, aset=aset, $
     fibermask=fibermask, bestcorr=bestcorr, $
     color=color
   if (NOT keyword_set(wset)) then return, 0
   traceset2xy, wset, xx, yy

   ;----------
   ; Compute fiber-to-fiber flat-field variations, and append to the end of
   ; the reduced flat file (HDU #5).  First see if we've already done this.

   fflat = mrdfits(flatname,5)
   if (NOT keyword_set(fflat)) then begin
      flat_flux = mrdfits(flatname,0)
      flat_ivar = mrdfits(flatname,1)
      fflat = fiberflat(flat_flux, flat_ivar, wset, fibermask=fibermask, $
       /dospline)
      flat_flux = 0 ; clear memory
      flat_ivar = 0 ; clear memory
      mwrfits, fflat, flatname
      fflat = 0 ; clear memory
   endif

   ;----------
   ; Write out wavelength solution

   mwrfits, wset, outarc ;, /create

   wavemin = 10^(min(yy))
   wavemax = 10^(max(yy))
   nlamps = (size(xpeak,/dimens))[1]
   rstruct = create_struct('FLAVOR', 'arc', $
                           'CAMERA', camname, $
                           'WAVEMIN', wavemin, $
                           'WAVEMAX', wavemax, $
                           'BESTCORR', bestcorr, $
                           'NLAMPS', nlamps )

   return, rstruct
end

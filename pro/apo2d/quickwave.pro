function quickwave, arcname, tsetfile, wsetfile, fflatfile, radius=radius

   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(wsetfile) NE 1) then return, 0
   if (n_elements(tsetfile) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0

   ;----------
   ; Read in image

   sdssproc, arcname, arcimg, hdr=hdr, color=color, camname=camname

   ;----------
   ; Read in the reduced flat

   tset = mrdfits(tsetfile,2)
   fibermask = mrdfits(tsetfile,4)
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
   ; Compute fiber-to-fiber flat-field variations.
   ; First see if we've already done this.

   fflatexist = keyword_set(findfile(fflatfile))
   if (NOT fflatexist) then begin
      flat_flux = mrdfits(tsetfile,0)
      flat_ivar = mrdfits(tsetfile,1)
      fflat = fiberflat(flat_flux, flat_ivar, wset, fibermask=fibermask, $
       /dospline)
      flat_flux = 0 ; clear memory
      flat_ivar = 0 ; clear memory
      mwrfits, fflat, fflatfile, /create
      fflat = 0 ; clear memory
   endif

   ;----------
   ; Write out wavelength solution

   mwrfits, wset, wsetfile, /create

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

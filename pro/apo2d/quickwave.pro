function quickwave, arcname, flatname, outarc, radius=radius

   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(outarc) NE 1) then return, 0
   if (n_elements(flatname) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0

   ; Read in image
   sdssproc, arcname, arcimg, hdr=hdr, color=color
   camname = strtrim(sxpar(hdr, 'CAMERAS'),2)

   tset = mrdfits(flatname,1)
   fibermask = mrdfits(flatname,3)

   traceset2xy, tset, ycen, xcen 

   ; Boxcar extract
   flux = extract_boxcar(arcimg, xcen, radius=radius)

   ; Estimate inverse variance
   fluxivar = 1.0 / (abs(flux) + 10.0)

   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, aset=aset, $
     fibermask=fibermask, bestcorr=bestcorr, $
     color=color

   ; Only output file if a wavelength solution was found
   if (keyword_set(wset)) then begin
      mwrfits, wset, outarc ;, /create

      traceset2xy, wset, xx, yy
      wavemin = 10^(min(yy))
      wavemax = 10^(max(yy))
      nlamps = (size(xpeak,/dimens))[1]
      rstruct = create_struct('FLAVOR', 'arc', $
                              'CAMERA', camname, $
                              'WAVEMIN', wavemin, $
                              'WAVEMAX', wavemax, $
                              'BESTCORR', bestcorr, $
                              'NLAMPS', nlamps )
   endif else begin
      rstruct = 0
   endelse

   return, rstruct
end

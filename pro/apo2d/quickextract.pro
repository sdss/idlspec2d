function quickextract, flatname, arcname, sciname, outname, $
 radius=radius, filtsz=filtsz

   if (n_params() LT 4) then begin
      print, 'Syntax - quickextract, flatname, arcname, sciname, outname'
      return, 0
   endif

   if (n_elements(flatname) NE 1) then return, 0
   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(sciname) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0
   if (NOT keyword_set(filtsz)) then filtsz = 25

   ;----------
   ; Read in the raw science image

   sdssproc, sciname, image, hdr=hdr, color=color
   camname = strtrim(sxpar(hdr, 'CAMERAS'),2)
   exptime = sxpar(hdr, 'EXPTIME')

   ;----------
   ; Read in the reduced data from the flat and arc

   tset = mrdfits(flatname,1)
   plugsort = mrdfits(flatname,2)
   fibermask = mrdfits(flatname,3)
   wset = mrdfits(arcname,1)

   traceset2xy, tset, ytemp, xcen
   traceset2xy, wset, ytemp, logwave

   ;----------
   ; Extract and sky-subtract the science image

   ; Boxcar extract - no scattered light correction!
   flux = extract_boxcar(image, xcen, radius=radius)

   ; Estimate fluxivar
   fluxivar = 1.0 / (abs(flux) + 10.0)

   ; Sky-subtract
   skystruct = skysubtract(flux, fluxivar, plugsort, wset, $
    objsub, objsubivar, iskies=iskies, fibermask=fibermask)

   ;----------
   ; Analyze spectra for the sky level and signal-to-noise

   ; Select wavelength range to analyze
   if (strmid(color,0,1) EQ 'b') then begin
      icolor = 1
      wrange = [4000,5500]
      snmag = 19.2
   endif else begin
      icolor = 3
      wrange = [6910,8500]
      snmag = 18.9
   endelse

   ; Find which fibers are sky fibers + object fibers
;   iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
;    AND plugsort.fiberid GT 0)
   iobj = where(strtrim(plugsort.objtype,2) NE 'SKY' $
    AND plugsort.fiberid GT 0)

   ; Compute average (but median-filtered) flux and signal-to-noise
   nfiber = (size(flux,/dimens))[1]
   meanflux = fltarr(nfiber)
   meansn = fltarr(nfiber)
   for ifib=0, nfiber-1 do begin
      ; Select unmasked pixels in the wavelength range for this fiber
      iwave = where(logwave[*,ifib] GT alog10(wrange[0]) $
                AND logwave[*,ifib] LT alog10(wrange[1]) $
                AND objsubivar[*,ifib] GT 0, nwave)
      if (nwave GT 0) then begin
         meanfluxvec = djs_median( flux[iwave,ifib], width=filtsz<nwave, $
          boundary='reflect' )
         meansnvec = djs_median( objsub[iwave,ifib] $
          * sqrt(objsubivar[iwave,ifib]), width=filtsz<nwave, $
          boundary='reflect' )
         meanflux[ifib] = djs_mean(meanfluxvec)
         meansn[ifib] = djs_mean(meansnvec)
      endif
   endfor

   ; Compute median of the sky fibers
   if (iskies[0] NE -1) then begin
      skylevel = median( meanflux[iskies] )
      if (exptime GT 0) then skylevel = skylevel / exptime
   endif else begin
      skylevel = 0
   endelse

   if (iobj[0] NE -1) then begin
;      snoise2 =  djs_mean( meansn[iobj] )^2
      coeffs = fitsn(plugsort[iobj].mag[icolor], meansn[iobj])
      snoise2 = 10^(2.0 * poly(snmag, coeffs)) ; The 2.0 is to square the S/N
   endif else begin
      snoise2 = 0
   endelse

   rstruct = create_struct('FLAVOR', 'science', $
                           'CAMERA', camname, $
                           'SKYPERSEC', skylevel, $
;                           'SN2VECTOR', meansn, $
                           'SN2', snoise2 )

   ;----------
   ; Write out the extracted spectra

   mwrfits, objsub, outname, /create
   mwrfits, ojbsubivar, outname
   mwrfits, meansn, outname

   return, rstruct
end


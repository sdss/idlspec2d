function quickextract, flatname, arcname, sciname, outname, radius=radius

   if (N_PARAMS() LT 4) then begin
      print, 'Syntax - quickextract, flatname, arcname, sciname, outname'
      return, 0
   endif

   if (n_elements(flatname) NE 1) then return, 0
   if (n_elements(arcname) NE 1) then return, 0
   if (n_elements(sciname) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0

   ; Read in image
   sdssproc, sciname, image, hdr=hdr, color=color
   camname = strtrim(sxpar(hdr, 'CAMERAS'),2)

   tset = mrdfits(flatname,1)
   plugsort = mrdfits(flatname,2)
   fibermask = mrdfits(flatname,3)
   wset = mrdfits(arcname,1)

   traceset2xy, tset, ytemp, xcen
   traceset2xy, wset, ytemp, logwave

   ; Boxcar extract
   flux = extract_boxcar(image, xcen, radius=radius)

   ; Estimate fluxivar
   fluxivar = 1.0/(abs(flux) + 10.0)

   mwrfits, flux, outname, /create
   mwrfits, fluxivar, outname

   ; Select wavelength range to analyze
   color = strmid(camname,0,1)
   if (color EQ 'b') then wrange = [4000,5500] $
    else wrange = [6910,8500]

   ; Find which fibers are sky fibers + object fibers
   isky = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
    AND plugsort.fiberid GT 0)
   iobj = where(strtrim(plugsort.objtype,2) NE 'SKY' $
    AND plugsort.fiberid GT 0)

   nfiber = (size(flux,/dimens))[1]
   medflux = fltarr(nfiber)
   medivar = fltarr(nfiber)
   for ifib=0, nfiber-1 do begin
      iwave = where(logwave[*,ifib] GT alog10(wrange[0]) $
                AND logwave[*,ifib] LT alog10(wrange[1]) )
      medflux[ifib] = median( flux[iwave,ifib] )
      medivar[ifib] = median( fluxivar[iwave,ifib] )
   endfor

   if (isky[0] NE -1) then begin
      tmpw = logwave[*,isky]
      tmpf = flux[*,isky]
      skylevel = djs_mean( medflux[isky] )
   endif else begin
      skylevel = 0
   endelse

   if (iobj[0] NE -1) then begin
      snoise =  djs_mean( (medflux[iobj]) * sqrt(medivar[iobj]) )
   endif else begin
      snoise = 0
   endelse

   rstruct = create_struct('FLAVOR', 'science', $
                           'CAMERA', camname, $
                           'SKYLEVEL', skylevel, $
                           'SN2', snoise^2 )

   return, rstruct
end


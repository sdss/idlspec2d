;+
; NAME:
;   quickextract
;
; PURPOSE:
;   Extraction with quick scattered light and sky subtraction.
;   S/N is estimated and output for html generation.
;
; CALLING SEQUENCE:
;   rstruct = quickextract(tsetfile, wsetfile, fflatfile, rawfile, outsci, $
;           [radius=radius, filtsz=filtsz])
;
; INPUTS:
;   tsetfile   - Name of fits file which contains matched trace
;   wsetfile   - Name of fits file which contains matched wavelengths
;   fflatfile  - Name of fits file which containes flat field vectors
;   rawfile    - Name of SDSS science frame to extract
;   outsci     - Name of fits file to store workings of quickextract
;
; OPTIONAL INPUTS:
;   radius     - Radius for boxcar extraction (default 3)
;   filtsz     - Median filter size to apply before average S/N is calculated.
;
; OUTPUT:
;   rstruct    - Results to be added html file upon completion
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   calcscatimage
;   divideflat
;   extract_image
;   extract_boxcar()
;   fileandpath()
;   find_whopping()
;   quickboxcar()
;   mrdfits()
;   mwrfits
;   sdssproc
;   skysubtract()
;   traceset2xy
;
; REVISION HISTORY:
;   3-Apr-2000  Written by S. Burles & D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function quickextract, tsetfile, wsetfile, fflatfile, rawfile, outsci, $
 radius=radius, filtsz=filtsz

   if (n_params() LT 4) then begin
      print, 'Syntax - quickextract, tsetfile, wsetfile, fflatfile, rawfile, outsci'
      return, 0
   endif

   if (n_elements(tsetfile) NE 1) then return, 0
   if (n_elements(wsetfile) NE 1) then return, 0
   if (n_elements(rawfile) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0
   if (NOT keyword_set(filtsz)) then filtsz = 25

   ;----------
   ; Read in the raw science image

   sdssproc, rawfile, image, invvar, hdr=hdr, color=color, camname=camname
   exptime = sxpar(hdr, 'EXPTIME')

   ;----------
   ; Read in the reduced data from the flat and arc

   tset = mrdfits(tsetfile,2)
   plugsort = mrdfits(tsetfile,3)
   fibermask = mrdfits(tsetfile,4)
   fflat = mrdfits(fflatfile,0)
   fflatmask = mrdfits(fflatfile,1)
   fibermask = fibermask OR fflatmask
   wset = mrdfits(wsetfile,1)
   traceset2xy, wset, ytemp, logwave

   ;---------------------------------------------------------------------------
   ; Extract and sky-subtract the science image
   ;---------------------------------------------------------------------------

   ;----------
   ; First let's try scattered light fit

   nrow = (size(image))[2]
   ncol = (size(image))[1]
   skiprow = 8
   yrow = lindgen(nrow/skiprow) * skiprow + skiprow/2
   nfirst = n_elements(yrow)
   wfixed = [1,1] ; Fit gaussian height + width

   sigma = 1.0
   proftype = 1 ; Gaussian
   npoly=8
   nterms=2

   traceset2xy, tset, ytemp, xcen

   ;----------
   ; Calculate the shift of the traces between the flat and science exposures

   xnew = match_trace(image, invvar, xcen)
   bestlag = median(xnew-xcen)
   minlag = min(xnew-xcen)
   maxlag = max(xnew-xcen)
   splog, 'Match_trace range: ', minlag, bestlag, maxlag
   if (bestlag LT -0.15 OR bestlag GT 0.15) then $
    splog, 'WARNING: Large flexure flat<->science (Post-calibs recommended!)'

   extract_image, image, invvar, xcen, sigma, tempflux, tempfluxivar, $
    proftype=proftype, wfixed=wfixed, yrow=yrow, highrej=5, lowrej=5, $
    npoly=npoly, ansimage=ansimage, relative=1

   ntrace = (size(tempflux))[2]
   scatfit = calcscatimage(ansimage[ntrace*nterms:*,*], yrow)
   scatflux = extract_boxcar(scatfit, xcen, radius=radius)

   exptime_factor = (exptime/900.0) > 1.0
   if (strmid(color,0,1) EQ 'b') then scatlimit = 20 *exptime_factor $
    else scatlimit = 30 * exptime_factor

   scatmed = median(scatfit)
   scatmax = max(scatfit)

   if (scatmed GT scatlimit) then $
     splog, 'WARNING: Scattered light median = ', scatmed, ' electrons' $
      + ' (WARM CCD OR TWILIGHT?)' $
    else $
     splog, 'Scattered light median = ', scatmed, ' electrons'
   if (scatmax GT 2*scatlimit) then $
     splog, 'WARNING: Scattered light max = ', scatmax, ' electrons' $
      + ' (WARM CCD OR TWILIGHT?)' $
    else $
     splog, 'Scattered light max = ', scatmax, ' electrons'

   ;----------
   ; Boxcar extract - no scattered light correction!

   flux = quickboxcar(image, invvar, tset=tset, fluxivar=fluxivar)
   fluxsub = flux - scatflux
   nfiber = (size(flux,/dimens))[1]

   ; Flat-field
   divideflat, fluxsub, invvar=fluxivar, fflat, /quiet

   ;--------------------------------------------------------------------
   ; Check for whopping fibers in SOS reductions, 
   ; especially useful for flagging affected sky fibers

   scrunch = djs_median(flux, 1) ; Find median counts/row in each fiber
   whopping = find_whopping(scrunch, 10000.0, whopct)
   if (whopping[0] NE -1) then begin
     splog, 'WARNING: Whopping fiber(s) at ', whopping, $
            '  (may have adverse affect on S/N)'
     wp = [whopping - 2 , whopping -1, whopping, whopping+1 , whopping+2]
     wp = (wp > 0) < (nfiber - 1)
     fibermask[wp] = fibermask[wp] OR pixelmask_bits('NEARWHOPPER')
   endif
   
   iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
      AND (plugsort.fiberid GT 0) AND (fibermask EQ 0), nskies)

   if nskies GT 10 then begin
      skylevel = djs_median(fluxsub[*,iskies], 1)
      outlier = (sort(skylevel))[[0,1,nskies-2,nskies-1]]
      fibermask[iskies[outlier]] = fibermask[iskies[outlier]] $
                                OR fibermask_bits('BRIGHTSKY')
   endif


   ; Sky-subtract
   skystruct = skysubtract(fluxsub, fluxivar, plugsort, wset, $
    objsub, objsubivar, iskies=iskies, fibermask=fibermask)

   ;---------------------------------------------------------------------------
   ; Analyze spectra for the sky level and signal-to-noise
   ;---------------------------------------------------------------------------

   ;----------
   ; Select wavelength range to analyze

   if (strmid(color,0,1) EQ 'b') then begin
      icolor = 1
      wrange = [4000,5500] ; coverage of g-band
      fitmag = [18.2, 19.7]
      snmag = 20.2
   endif else begin
      icolor = 3
      wrange = [6910,8500] ; coverage of i-band
      fitmag = [17.9, 19.4]
      snmag = 19.9
   endelse

   ;----------
   ; Find which fibers are sky fibers + object fibers

   iobj = where(strtrim(plugsort.objtype,2) NE 'SKY' $
    AND plugsort.fiberid GT 0)

   ;----------
   ; Compute average (but median-filtered) flux and signal-to-noise

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

   ;----------
   ; Compute median of the sky fibers

   if (iskies[0] NE -1) then begin
      skylevel = median( meanflux[iskies] )
      if (exptime GT 0) then skylevel = skylevel / exptime
   endif else begin
      skylevel = 0
   endelse

   if (iobj[0] NE -1) then begin
;      fitmag = snmag+[-1.0,0.5]
      if (n_elements(where(plugsort[iobj].mag[icolor] GT fitmag[0] $
       AND plugsort[iobj].mag[icolor] LT fitmag[1])) LT 20) $
       then fitmag[0] = 1.0

      coeffs = fitsn(plugsort[iobj].mag[icolor], meansn[iobj], fitmag=fitmag)
      if (keyword_set(coeffs)) then $
       snoise2 = 10^(2.0 * poly(snmag, coeffs)) $ ; The 2.0 is to square the S/N
      else $
       snoise2 = 0
   endif else begin
      snoise2 = 0
   endelse

   rstruct = create_struct('SCIFILE', fileandpath(outsci), $
                           'SKYPERSEC', skylevel, $
;                           'MAG', plugsort.mag[icolor], $
                           'SN2VECTOR', meansn^2, $
                           'SN2', snoise2 )

   ;----------
   ; Write out the extracted spectra

   mwrfits, objsub, outsci, /create
   mwrfits, objsubivar, outsci
   mwrfits, meansn, outsci

   return, rstruct
end


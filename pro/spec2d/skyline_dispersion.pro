;+
; NAME:
;   skyline_dispersion
;
; PURPOSE:
;   Measure and apply broadening to dispersion solution from an arc traceset
;   to agree with widths of sky lines, which may be broader due to flexure.
;
; CALLING SEQUENCE:
;   skydispset = skyline_dispersion(flux, fluxivar, xcen, iskies, dispset)
;
; INPUTS:
;   flux        - Object extracted flux [NPIX,NTRACE]
;   fluxivar    - Corresponding inverse variance [NPIX,NTRACE]
;   xcen        - xpeaks of sky lines [NTRACE,NLINES]
;   iskies      - Fiber indices corresponding to good sky fibers;
;                 default to using all fibers if this is not set.
;   dispset     - Original arcline solution for resolution vs. pixel
;                 (typically measured in native pixels)
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   skydispset  - Modified version of DISPSET with dispersions broadened
;                 to agree with the sky line widths.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The pixel dispersion can be broadened, but never made smaller.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   extract_image
;   splog
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   20-Feb-2002  Written by S. Burles, MIT (Adapted from fitdispersion)
;-
;------------------------------------------------------------------------------
function skyline_dispersion, flux, fluxivar, xcen, iskies, dispset

   ndim = size(xcen, /n_dimen)
   dims = size(xcen, /dimens)
   if (ndim EQ 1) then begin
      ntrace = 1
      nline = dims[0]
   endif else begin
      ntrace = dims[0]
      nline = dims[1]
   endelse
   if (n_elements(iskies) EQ 0) then iskies = lindgen(ntrace)
   npix = (size(flux,/dimen))[0]
   xmin = 0.0
   xmax = npix - 1.0
   nskies = n_elements(iskies)
   splog, 'Number of sky lines = ', nline
   splog, 'Number of fibers = ', ntrace
   splog, 'Number of sky fibers = ', nskies

   ;----------
   ; Sort sky-line positions, which must be done for EXTRACT_IMAGE.

   isort = sort(xcen[0,*])
   xsky = xcen[*,isort]
   ysky = rebin(lindgen(ntrace), ntrace, nline)

   ;----------
   ; Construct a mask that is nonzero only for those pixels within
   ; +/- 12 pixels from each sky line center on each fiber.

   smask = make_array(size=size(flux), /byte)
   for itrace=0, ntrace-1 do begin
      for iline=0, nline-1 do begin
         i1 = floor(xsky[itrace,iline] - 12) > 0
         i2 = ceil(xsky[itrace,iline] + 12) < (npix - 1)
         if (i1 LT npix-1 AND i2 GT 0) then smask[i1:i2,itrace] = 1
      endfor
   endfor

   ;----------
   ; Measure the arc-line dispersions at the location of each sky line

   traceset2xy, dispset, transpose(xsky), arcwidth
   sigma = transpose(arcwidth)

   ;----------
   ; Extract sky lines from the [NPIX,NTRACE] image, measuring the
   ; wavelength sigmas during extraction for every fiber (sky and non-sky).
   ; SKYLINEFLUX = [NFIBER,NLINE]
   ; SKYLINEIVAR = [NFIBER,NLINE]

   extract_image, flux, fluxivar*smask, xsky, transpose(arcwidth), $
    skylineflux, skylineivar, ansimage=ansimage, wfixed=[1,1], $
    highrej=10, lowrej=10, relative=1, npoly=5, proftype=1

   ;----------
   ; Compare the width terms for arcs and sky lines in only the sky fibers
   ; GMASK    = [NLINE,NFIBER]
   ; SKYGMASK = [NLINE,NFIBER] --> Same as above, but nonzero only for skies
   ; ARCWIDTH = [NLINE,NFIBER]
   ; SKYWIDTH = [NLINE,NFIBER]

   gmask = transpose(skylineivar) GT 0 AND transpose(skylineflux) GT 0

   skywidth = arcwidth * (1 + ansimage[lindgen(nline)*2+1,*] $
    / (transpose(skylineflux) * gmask + (1-gmask)) )
   quad_diff = skywidth^2 - arcwidth^2

   skygmask = 0 * gmask
   skygmask[*,iskies] = gmask[*,iskies]

   ;----------
   ; Log the differences between the arc and sky-line widths for each line
   ; (only using sky fibers).

   for iline=0, nline-1 do begin
      igood = where(skygmask[iline,iskies] GT 0, ngood)
      if (ngood GT 0) then begin
         djs_iterstat, arcwidth[iline,igood], median=med1, sigma=sig1
         djs_iterstat, skywidth[iline,igood], median=med2, sigma=sig2
         djs_iterstat, quad_diff[iline,igood], median=med3, sigma=sig3
         splog, 'Line #', iline, ' median sigma(arc)=', med1, $
          ' sigma(sky)=', med2
         splog, 'Line #', iline, ' sigma^2(add)=', med3, ' +/- ', sig3
      endif
   endfor

   ;----------
   ; Measure the median value of the quadrature difference between
   ; the sky-line widths and arc-line widths (only using sky fibers).

   igood = where(skygmask GT 0, ngood)
   if (ngood GT 0) then begin
      djs_iterstat, quad_diff[igood], median=addsig2, sigma=addsig2_err
      splog, 'Flexure sigma^2(add) = ', addsig2, ' +/- ', addsig2_err
      if (addsig2 LT 0) then $
       splog, 'Warning: Truncating sigma to add in quadrature to zero'
   endif else begin
      splog, 'Flexure sigma^2(add) = (undefined)'
   endelse

   ;----------
   ; Re-fit the modified line-dispersion vectors, using the same
   ; number of coefficients as for the input trace set.

   if (keyword_set(addsig2)) then begin
      traceset2xy, dispset, xtmp, newsigma
      newsigma = sqrt(newsigma^2 + (addsig2 > 0))

      ncoeff = (size(dispset.coeff, /dimens))[0]
      xy2traceset, xtmp, newsigma, skydispset, ncoeff=ncoeff, $
       xmin=dispset.xmin, xmax=dispset.xmax, yfit=yfit
   endif else begin
      skydispset = dispset
   endelse

   return, skydispset
end
;------------------------------------------------------------------------------

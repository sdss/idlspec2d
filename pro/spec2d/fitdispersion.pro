;+
; NAME:
;   fitdispersion
;
; PURPOSE:
;   Fit polynomials to the line width (dispersion axis) for each fiber bundle
;
; CALLING SEQUENCE:
;   dispset = fitdispersion(arc_flux, arc_fluxivar, xcen, $
;    [ ncoeff=, sigma=, xmin=, xmax=, medwidth= ] )
;
; INPUTS:
;   arc_flux     - arc image extracted flux
;   arc_fluxivar - corresponding inverse variance
;   xcen         - xpeaks of arc lines [ntrace,nlines]
;
; OPTIONAL KEYWORDS:
;   ncoeff     - Order of legendre polynomial to apply to width vs. row;
;                default to 4.
;   sigma      - The SIGMA input to EXTRACT_IMAGE when determining ANSIMAGE;
;                default to 1.0 pix.  This can be a scalar, an [NFIBER] vector,
;                or an [NLINE,NFIBER] array.
;   xmin       - Lowest row number for trace set; default to 0.
;   xmax       - Highest row number for trace set; default to 2047.
;
; OUTPUTS:
;   dispset   - Traceset structure containing fit coefficients
;
; OPTIONAL OUTPUTS:
;   medwidth  - Median dispersion widths in each of the 4 quadrants
;               of the CCD, ordered LL,LR,UL,UR.
;
; COMMENTS:
;   Used to fill arcstruct.dispset, which can then be applied
;   to PSF-corrected sky subtraction.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;   xy2traceset
;
; REVISION HISTORY:
;   01-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitdispersion, arc_flux, arc_fluxivar, xcen_inp, $
 sigma=sigma, ncoeff=ncoeff, xmin=xmin, xmax=xmax, medwidth=medwidth

   if (NOT keyword_set(sigma)) then sigma = 1.0
   if (NOT keyword_set(ncoeff)) then ncoeff = 4
   if (NOT keyword_set(xmin)) then xmin = 0.0
   if (NOT keyword_set(xmax)) then xmax = 2047.0

   nline = (size(xcen_inp,/dimen))[1]
   ntrace = (size(xcen_inp,/dimen))[0]
   if (ntrace NE 320) then $
    message, 'Must have 320 traces!'
   npix = (size(arc_flux,/dimen))[0]

   if (nline LT 3) then begin
      splog, 'WARNING: Too few lines for dispersion traceset (', nline, ')'
      return, -1
   endif

   ;----------
   ; Sort XCEN, which is necessary for the call to EXTRACT_IMAGE.

   isort = sort(xcen_inp[0,*])
   xcen = xcen_inp[*,isort]

   ;----------
   ; Construct a mask that is nonzero only for those pixels within
   ; +/- 12 pixels from each arc line center on each fiber.

   arcmask = make_array(size=size(arc_flux), /byte)
   for itrace=0, ntrace-1 do begin
      for iline=0, nline-1 do begin
         i1 = floor(xcen[itrace,iline] - 12) > 0
         i2 = ceil(xcen[itrace,iline] + 12) < (npix - 1)
         if (i1 LT npix-1 AND i2 GT 0) then arcmask[i1:i2,itrace] = 1
      endfor
   endfor

   ;----------
   ; Extract arc lines from the [NPIX,NTRACE] image, measuring the
   ; wavelength sigmas during extraction for every fiber.

   extract_image, arc_flux, arc_fluxivar*arcmask, xcen, sigma, $
    arclineflux, arclineivar, ansimage=ansimage, wfixed=[1,1], $
    highrej=10, lowrej=10, relative=1, npoly=5

   ;----------
   ; Determine the widths from the output array from EXTRACT_IMAGE.

   gmask = transpose((arclineivar GT 0) * (arclineflux GT 0))
   igood = where(gmask)
   widthterm = ansimage[lindgen(nline)*2+1,*]
   width = make_array(size=size(transpose(arclineflux)), /float)
   if (igood[0] NE -1) then $
    width[igood] = (1 + widthterm[igood] / (transpose(arclineflux))[igood])

   ndim = size(sigma, /n_dimen)
   if (n_elements(sigma) EQ 1) then begin
      width = width * sigma[0]
   endif else if (ndim EQ 1) then begin
      for itrace=0, ntrace-1 do $
       width[*,itrace] = width[*,itrace] * sigma[itrace]
   endif else if (ndim EQ 2) then begin
      if (n_elements(sigma) NE n_elements(width)) then $
       message, 'Dimensions of SIGMA and WIDTH do not agree'
      width = width * sigma
   endif else begin
      message, 'Unsupported number of elements for SIGMA'
   endelse

   ;----------
   ; Perform median across bundles on good arclines only
   ; somewhat tedious, but it works

   width = reform(width,nline,20,16)
   gmask = reform(gmask,nline,20,16)
   width_bundle = fltarr(nline,16)

   for iline=0, nline-1 do begin
     for j=0, 15 do begin
        ss = where(gmask[iline,*,j] AND width[iline,*,j] GT 0)
        if (ss[0] NE -1) then $
         width_bundle[iline,j] = djs_median(width[iline,ss,j]) 
     endfor
   endfor

   width_final = rebin(width_bundle, nline, ntrace, /sample)

   ;----------
   ; Turn the widths back into a traceset.

   xy2traceset, transpose(xcen), width_final, dispset, $
    ncoeff=ncoeff, xmin=xmin, xmax=xmax

   ;----------
   ; Trigger warning messages if widths are too large.

   ; Which lines are on the bottom/top of the CCD?
;   ilo = where(xcen[0,*] LE 0.5*xmax) > 0
;   ihi = where(xcen[0,*] GT 0.5*xmax) > 0
   traceset2xy, dispset, xx, width_fit
   medwidth = [ median(width_fit[0:npix/2-1,0:ntrace/2-1]), $
                median(width_fit[0:npix/2-1,ntrace/2:ntrace-1]), $
                median(width_fit[npix/2:npix-1,0:ntrace/2-1]), $
                median(width_fit[npix/2:npix-1,ntrace/2:ntrace-1]) ]

   splog, ((max(medwidth) GT 1.10) ? 'WARNING: ' : '') $
    + 'Median wavelength widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (LL LR UL UR)'

   return, dispset
end
;------------------------------------------------------------------------------

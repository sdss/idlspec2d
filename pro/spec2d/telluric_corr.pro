;+
; NAME:
;   telluric_corr
;
; PURPOSE:
;   Use SPECTROPHOTO and REDDEN_STD's to fit telluric features.
;
; CALLING SEQUENCE:
;   telluric_factor = telluric_corr(flux, fluxivar, wset, plugsort, $
;    fibermask=, tellbands=, ncoeff=, pixspace=, /dospline, $
;    nord=, upper=, lower=, plottitle=
;
; INPUTS:
;   flux         - Sky-subtracted extracted spectra [NX,NTRACE]
;   fluxivar     - Inverse variance for FLUX [NX,NTRACE]
;   wset         - Wavelength coefficients as a trace set
;   plugsort     - Plugmap entries
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   tellbands  - Structure array defining the telluric bands and the
;                wavelengths used to fit the continuum.  If not set,
;                then default values are used.  See the code for the
;                definition of this structure.
;   ncoeff     - Number of coefficients used in constructing continuum;
;                default to 4.
;   pixspace   - Approximate spacing in pixels for break points in the
;                spline fits to individual fibers; default to 10 pixels.
;   dospline   - If this keyword is set, then fit the continuum
;                to splines (using PIXSPACE) rather than to a Legendre
;                polynomial (using NCOEFF).
;   plottitle  - Title for QA plot; if not set, then do not plot.
;
; PARAMETERS FOR SLATEC_SPLINEFIT:
;   nord
;   upper
;   lower
;
; OUTPUTS:
;   telluric_factor - Telluric correction for each pixel in flux array;
;                     set to 1 where there is no correction.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Objects that are SPECTROPHOTO_STD or REDDEN_STD (as listed in the plugmap
;   structure) are selected for constructing the telluric-correction factor.
;
; EXAMPLES:
;
; BUGS:
;   Other fainter telluric bands at [7178,7212], [8130,8206].
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   djs_reject()
;   fibermask_bits()
;   func_fit()
;   bspline_valu()
;   bspline_iterfit()
;   splog
;   traceset2xy
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;   15-Aug-2000  Major modifications by D. Schlegel to fit all telluric
;                bands with one spline, and to do this fit on the spectra
;                before flux-correction.
;-
;------------------------------------------------------------------------------
function telluric_corr,flux, fluxivar, wset, plugsort, $
 fibermask=fibermask, tellbands=tellbands, ncoeff=ncoeff, pixspace=pixspace, $
 nord=nord, upper=upper, lower=lower, dospline=dospline, $
 plottitle=plottitle

   ndim = size(flux, /n_dimen)
   if (ndim NE 2) then message, 'Expecting 2-D flux array'
   dims = size(flux, /dimens)
   npix = dims[0]
   ntrace = dims[1]

   if (N_elements(pixspace) EQ 0) then pixspace = 10
   if (N_elements(ncoeff) EQ 0) then ncoeff = 4
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(ntrace)

   tellcorr = fltarr(npix,ntrace) + 1.0 ; Return value if this routine fails

   ;----------
   ; Construct a structure describing which pixels to use for continuum-fitting,
   ; and which to use to construct the telluric-correction.

   if (NOT keyword_set(tellbands)) then begin
      tellbands1 = { TELLBAND, $
       twave1: 6850., $
       twave2: 6960., $
       cwave1: [6600., 6950., 0], $
       cwave2: [6860., 7200., 0] }
      tellbands2 = { TELLBAND, $
       twave1: 7150., $
       twave2: 7350., $
       cwave1: [7050., 7115., 7340.], $
       cwave2: [7160., 7130., 7440.] }
      tellbands3 = { TELLBAND, $
       twave1: 7560., $
       twave2: 7720., $
       cwave1: [7400., 7700., 0], $
       cwave2: [7580., 8000., 0] }
      tellbands4 = { TELLBAND, $
       twave1: 8105., $
       twave2: 8240., $
       cwave1: [8000., 8225., 0], $
       cwave2: [8105., 8325., 0] }
;      tellbands5 = { TELLBAND, $
;       twave1: 8530., $
;       twave2: 8865., $
;       cwave1: [8490., 8865., 0], $
;       cwave2: [8530., 8905., 0] }
;      tellbands6 = { TELLBAND, $
;       twave1: 8644., $
;       twave2: 8697., $
;       cwave1: [8604., 8697., 0], $
;       cwave2: [8644., 8737., 0] }
      tellbands = [tellbands1, tellbands2, tellbands3, tellbands4]
   endif

   ;----------
   ; Use SPECTROPHOTO_STD and REDDEN_STD to correct telluric absorption

   qphoto = strtrim(plugsort.objtype,2) EQ 'SPECTROPHOTO_STD' $
        OR strtrim(plugsort.objtype,2) EQ 'REDDEN_STD'
   qgfib = ((fibermask AND fibermask_bits('NOPLUG')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADTRACE')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADFLAT')) EQ 0) $
       AND ((fibermask AND fibermask_bits('NEARWHOPPER')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADARC')) EQ 0)
   tindx = where(qphoto AND qgfib, ntell)

   splog, 'Number of telluric standards = ', ntell

   if (ntell EQ 0) then begin
      print, 'WARNING: No telluric correction stars'
      return, tellcorr
   endif

   ;----------
   ; Compute the wavelengths for all fibers from the trace set.

   traceset2xy, wset, pixnorm, loglam

   ;----------

   if (keyword_set(dospline)) then begin
      ;----------
      ; Always select the same break points in log-wavelength for all fibers

      nbkpts = fix(npix / pixspace) + 2
      bkpt = findgen(nbkpts) * (max(loglam) - min(loglam)) / (nbkpts-1) $
       + min(loglam)
   endif

   ;----------
   ; Initially, set the telluric correction at all wavelengths to unity
   ; with a S/N of 100.
   ; In the telluric bands, these values will be replaced with the measured
   ; values.  Then all wavelengths will be fit with a spline.  This constrains
   ; the spline to go to unity outside of the telluric bands.

   fitflux = fltarr(npix,ntell) + 1
   fitivar = fltarr(npix,ntell) + 1.0e4

   ;---------------------------------------------------------------------------
   ; LOOP OVER EACH TELLURIC BAND
   ;---------------------------------------------------------------------------

   for iband=0, n_elements(tellbands)-1 do begin

      ;----------
      ; Set masks to 1 where continuum and telluric bands are defined

      cmask = bytarr(npix,ntrace)
      tmask = bytarr(npix,ntrace)

      for i=0, n_elements(tellbands[iband].cwave1)-1 do $
       cmask = cmask $
        OR (loglam GE alog10(tellbands[iband].cwave1[i]) $
        AND loglam LE alog10(tellbands[iband].cwave2[i]))
      for i=0, n_elements(tellbands[iband].twave1)-1 do $
       tmask = tmask $
        OR (loglam GE alog10(tellbands[iband].twave1[i]) $
        AND loglam LE alog10(tellbands[iband].twave2[i]))

      ;----------
      ; Create sub-arrays that contain only the telluric stars

      tellloglam = loglam[*,tindx]
      tellflux = flux[*,tindx]
      tellivar = fluxivar[*,tindx]
      tellcmask = cmask[*,tindx]
      telltmask = tmask[*,tindx]

      ;----------
      ; Fit continuum to each telluric standard

      nuse = 0

      for itell=0, ntell-1 do begin

         ;----------
         ; Insist that at least 70% of the pixels in both the continuum regions
         ; and the telluric regions are unmasked, and that there are at least
         ; 2 pixels of each.

         cnum = total(tellcmask[*,itell] * (tellivar[*,itell] NE 0))
         cfrac = cnum / (total(tellcmask[*,itell]) > 1)
         tnum = total(telltmask[*,itell] * (tellivar[*,itell] NE 0))
         tfrac = tnum / (total(telltmask[*,itell]) > 1)

         indc = where(tellcmask[*,itell] NE 0)
         indt = where(telltmask[*,itell] NE 0)

         if (cnum GE 2 AND cfrac GE 0.70 $
          AND tnum GE 2 AND tfrac GE 0.70) then begin

            if (keyword_set(dospline)) then begin

               ;---------------------------------------------------------------
               ; SPLINE FIT TO CONTINUUM
               ;---------------------------------------------------------------

               istart = (where(bkpt GT min(tellloglam[indc,itell])))[0]
               istart = (istart - 1) > 0
               iend = (where(bkpt GT max(tellloglam[indc,itell])))[0]
               if (iend EQ -1) then iend = nbkpts-1

               ; Reject at most 10% of the points
               maxrej = ceil(0.10 * (istart - iend + 1))

               sset = bspline_iterfit(tellloglam[indc,itell], $
                tellflux[indc,itell], invvar=tellivar[indc,itell], $
                maxiter=10, upper=upper, lower=lower, $
                nord=nord, bkpt=bkpt[istart:iend], maxrej=maxrej, requiren=1)

               continuum = bspline_valu(tellloglam[*,itell], sset)

            endif else begin

               ;---------------------------------------------------------------
               ; LEGENDRE FIT TO CONTINUUM
               ;---------------------------------------------------------------

               ; Rescale X axis to be between 0 and 1 in the fitting region.
               xmin = min(tellloglam[indc,itell], max=xmax)
               xfit = (tellloglam[*,itell] - xmin) / (xmax - xmin)

               ; Iterate the fit using rejection at the 6-sigma level.
               ; At most, reject 25% of the points.
               iiter = 0
               maxiter = 5
               qdone = 0
               inmask = tellivar[indc,itell] NE 0
               outmask = 0
               while (NOT keyword_set(qdone) AND iiter LE maxiter) do begin
                  res = func_fit(xfit[indc], tellflux[indc,itell], $
                   ncoeff, invvar=tellivar[indc,itell], $
                   function_name='flegendre')
                  continuum = flegendre(xfit, ncoeff) # res

                  qdone = djs_reject(continuum[indc] - tellflux[indc,itell], $
                   tellivar[indc,itell], lower=6, upper=6, $
                   maxrej=ceil(0.25*n_elements(indc)), $
                   inmask=inmask, outmask=outmask)
                  iiter = iiter + 1
               endwhile

            endelse

            nuse = nuse + 1
            fitflux[indt,itell] = tellflux[indt,itell] / continuum[indt]
            fitivar[indt,itell] = tellivar[indt,itell] * (continuum[indt])^2
         endif else begin
            fitivar[*,itell] = 0
         endelse
      endfor

      splog, 'Use ', nuse, ' stars for telluric band ', $
       tellbands[iband].twave1, tellbands[iband].twave2
      if (nuse EQ 0) then $
       splog, 'WARNING: No stars for telluric band ', $
        tellbands[iband].twave1, tellbands[iband].twave2

   endfor

   ;----------
   ; Select the data points that describe the telluric factor, and
   ; sort by wavelength.

   indx = where(fitivar NE 0)
   isort = sort(tellloglam[indx])
   fitloglam = tellloglam[indx[isort]]
   fitflux = fitflux[indx[isort]]
   fitivar = fitivar[indx[isort]]

   ;----------
   ; Now bspline the features in contflux

   tellset = bspline_iterfit(fitloglam, fitflux, $
    nord=nord, maxiter=10, lower=lower, upper=upper, invvar=fitivar, $
    everyn=2*ntell, rejper=0.1, requiren=1)

   tellcorr = bspline_valu(loglam, tellset)

   ;---------------------------------------------------------------------------
   ; QA PLOT
   ;---------------------------------------------------------------------------

   if (keyword_set(plottitle)) then begin

      nbands = n_elements(tellbands)
      if (nbands LE 3) then $
       !p.multi = [0,1,nbands] $
      else $
       !p.multi = [0,2,ceil(nbands/2.)]

      xmargin = 50 ; Extra plot range in Angstroms

      for iband=0, n_elements(tellbands)-1 do begin

         xlo = min(tellbands[iband].twave1) - xmargin
         xhi = max(tellbands[iband].twave2) + xmargin
         indx = where(10^fitloglam GE xlo AND 10^fitloglam LE xhi)
         indx = indx[ sort(fitloglam[indx]) ] ; Sort the wavelengths
         if (indx[0] NE -1) then begin
            ; Plot the data points used for the fit
            djs_plot, 10^fitloglam[indx], fitflux[indx], ps=3, $
             xrange=[xlo,xhi], yrange=[0.0,1.5], ymargin=[2,4], $
             xstyle=1, xtitle='\lambda [A]', ytitle='Correction factor'

            if (iband EQ 0) then xyouts, 0.5, 0.95, plottitle, $
             charsize=1.5, align=0.5, /normal

            ; Overplot the fit
            tellplot = bspline_valu(fitloglam[indx], tellset)
            djs_oplot, 10^fitloglam[indx], tellplot
         endif

      endfor
      !p.multi = 0
   endif

   return, tellcorr
end
;------------------------------------------------------------------------------

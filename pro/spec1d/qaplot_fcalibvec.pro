;+
; NAME:
;   qaplot_fcalibvec
;
; PURPOSE:
;   Generate QA plot: apparent flux-calibration errors
;
; CALLING SEQUENCE:
;   qaplot_fcalibvec, loglam, objflux, objivar, synflux, plugmap, zans, $
;    [plottitle= ]
;
; INPUTS:
;   loglam     - Wavelengths in log-10(Angstroms) [NPIX]
;   objflux    - Object flux [NPIX,NOBJ]
;   objivar    - Object inverse variance [NPIX,NOBJ]
;   synflux    - Best-fit model flux [NPIX,NOBJ]
;   plugmap    - Plug-map structure [NOBJ]
;   zans       - Structure with redshift-fit information [NOBJ]
;
; OPTIONAL KEYWORDS:
;   plottitle  - TITLE of plot
;
; OUTPUTS:
;   Output plots only
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Plot the ratio of the observed spectra divided by the best-fit
;   model spectra as a function of observed wavelength.  Deviations
;   from unity are an empirical measure of the flux-calibration errors.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   djs_oplot
;   djs_xyouts
;
; REVISION HISTORY:
;   26-Feb-2002 Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function medfcalibvec, objflux, objivar, synflux, ifiber, minpts=minpts

   if (NOT keyword_set(minpts)) then minpts = 3
   npix = (size(objflux, /dimens))[0]
   medvec = fltarr(npix)

   if (ifiber[0] EQ -1) then return, medvec

   ;----------
   ; At each wavelength, select points where:
   ;   flux > 5 * flux_error
   ;   flux_synthetic > 0
   ;   flux_synthetic > 5 * flux_error

   qgood = objflux * sqrt(objivar) GT 5.0 AND synflux GT 0 $
    AND synflux * sqrt(objivar) GT 10.0

   for ipix=0, npix-1 do begin
      igood = where(objflux[ipix,ifiber] * sqrt(objivar[ipix,ifiber]) GT 5.0 $
       AND synflux[ipix,ifiber] GT 0 $
       AND synflux[ipix,ifiber] * sqrt(objivar[ipix,ifiber]) GT 5.0, ngood)
      if (ngood GE minpts) then medvec[ipix] = $
       median(objflux[ipix,ifiber[igood]] / synflux[ipix,ifiber[igood]])
   endfor

   ; Interpolate over bad values
   medvec = djs_maskinterp(medvec, medvec EQ 0, /const)

   ; Median-filter this vector, then smooth...
   medvec = median(medvec, 5)
   medvec = smooth(medvec, 5, /edge_truncate)

   return, medvec
end

;------------------------------------------------------------------------------
pro qaplot_fcalibvec, loglam, objflux, objivar, synflux, plugmap, zans, $
 plottitle=plottitle

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   if (ndim EQ 1) then begin
      nfiber = 1
      npix = dims[0]
   endif else begin
      nfiber = dims[0]
      npix = dims[1]
   endelse
   wave = 10^loglam
   xrange = minmax(wave) + [-100,100]
   csize = 1.5

   pmulti = !p.multi
   !p.multi = 0

   for specid=1, 2 do begin

      ; Set up plot
      plot, xrange, [1,1], xrange=xrange, yrange=[0.80,1.80], $
       /xstyle, /ystyle, xtitle='Wavelength [Ang]', xticklen=1.0, $
       ytitle='Observed/Model Flux + offset', title=plottitle, charsize=csize
      xpos = 0.90 * !x.crange[0] + 0.10 * !x.crange[1]
      xyouts, xpos, 1.75, 'Spectrograph #' + string(specid,format='(i1)'), $
       charsize=csize

      ; Plot residuals from QSOs
      ifiber = where(plugmap.spectrographid EQ specid $
       AND strtrim(zans.class,2) EQ 'QSO' $
       AND zans.zwarning EQ 0, nfiber)
      medvec = medfcalibvec(objflux, objivar, synflux, ifiber, minpts=1)
      djs_oplot, wave, medvec+0.20, color='blue'
      djs_oplot, !x.crange, [1.20,1.20]
      djs_xyouts, xpos, 1.65, color='blue', charsize=csize, $
       string(nfiber, format='("QSOs (", i3, ")")')

      ; Plot residuals from all fibers
      ifiber = where(plugmap.spectrographid EQ specid $
       AND zans.zwarning EQ 0, nfiber)
      medvec = medfcalibvec(objflux, objivar, synflux, ifiber)
      djs_oplot, wave, medvec+0.00
      djs_oplot, !x.crange, [1.00,1.00]
      djs_xyouts, xpos, 1.60, charsize=csize, $
       string(nfiber, format='("All good fibers (", i3, ")")')

      ;----------
      ; Trigger warning messages based upon the correction vector
      ; from all fibers.

      rmscorr = stddev(medvec,/double)
      outstring = string(specid, rmscorr, $
       format='("Spectro-",i1," RMS flux-corr errors = ", f6.3)')
      splog, outstring

      mincorr = min(medvec, max=maxcorr)
      outstring = string(specid, mincorr, maxcorr, $
       format='("Spectro-",i1," min/max flux-corr errors = ", f6.3, " ", f6.3)')
      splog, outstring
      if (mincorr LT 0.90 OR maxcorr GT 1.10) then $
       splog, 'WARNING: ' + outstring

      ; Plot residuals from F stars
      ifiber = where(plugmap.spectrographid EQ specid $
       AND (strtrim(plugmap.objtype,2) EQ 'SPECTROPHOTO_STD' $
           OR strtrim(plugmap.objtype,2) EQ 'REDDEN_STD') $
       AND zans.zwarning EQ 0, nfiber)
      medvec = medfcalibvec(objflux, objivar, synflux, ifiber)
      djs_oplot, wave, medvec+0.40, color='red'
      djs_oplot, !x.crange, [1.40,1.40]
      djs_xyouts, xpos, 1.70, color='red', charsize=csize, $
       string(nfiber, format='("Spectro-photo + reddening stars (", i3, ")")')
   endfor

   !p.multi = pmulti

   return
end
;------------------------------------------------------------------------------

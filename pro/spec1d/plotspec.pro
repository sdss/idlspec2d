;+
; NAME:
;   plotspec
;
; PURPOSE:
;   Plotting tool for a single fiber from Princeton Spectro-1D.
;
; CALLING SEQUENCE:
;   plotspec, plate, fiberid, [ mjd=, nsmooth= ]
;
; INPUTS:
;   plate      - Plate number
;   fiberid    - Fiber number
;
; OPTIONAL INPUTS:
;   mjd        - MJD number; if not set, then select the most recent
;                data for this plate (largest MJD).
;   nsmooth    - If set, then boxcar smooth both the object and synthetic
;                spectra with a width equal to NSMOOTH.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $SPECTRO_DATA/$PLATE/spPlate-$PLATE-$MJD.fits
;   $SPECTRO_DATA/$PLATE/spZbest-$PLATE-$MJD.fits
;   $IDLSPEC2D_DIR/templates/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   readspec
;   soplot
;   splot
;   sdss_flagname()
;   sxyouts
;   synthspec()
;   textoidl()
;
; REVISION HISTORY:
;   01-Sep-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro plotspec, plate, fiberid, mjd=mjd, nsmooth=nsmooth

   readspec, plate, fiberid, mjd=mjd, flux=objflux, flerr=objerr, $
    loglam=loglam, plug=plug, zans=zans
   wave = 10^loglam
   synflux = synthspec(zans, loglam=loglam)

   if (keyword_set(nsmooth)) then begin
      objflux = smooth(objflux, nsmooth)
      synflux = smooth(synflux, nsmooth)
   endif

   primtarget = sdss_flagname('TARGET', plug.primtarget, /concat)
   sectarget = sdss_flagname('TTARGET', plug.sectarget, /concat)

   csize = 2
   yrange = minmax(synflux)
   if (yrange[0] EQ yrange[1]) then yrange = minmax(objflux)
   ymin = 1.2 * yrange[0] - 0.2 * yrange[1] < 0
   ymax = -0.2 * yrange[0] + 1.2 * yrange[1]
   if (ymax EQ ymin) then ymax = ymin + 1

   title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(zans.mjd),2)
   splot, wave, objflux, yrange=[ymin,ymax], $
    xtitle='Wavelength [Ang]', $
    ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
    title=title, charsize=csize
   soplot, wave, objerr, color='red'
   soplot, wave, synflux, color='blue', lw=2

   xpos = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
   dypos = 0.05 * (!y.crange[0] - !y.crange[1])
   ypos = !y.crange[1] + 1.5 * dypos
   sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass $
    + '  z=' + strtrim(string(zans.z),2), $
    charsize=csize

   ypos = ypos + dypos
   sxyouts, xpos, ypos, $
    TeXtoIDL('X^2_r =' + strtrim(string(zans.rchi2), format='(f6.2)'),2)), $
    charsize=csize

   if (keyword_set(primtarget)) then begin
      ypos = ypos + dypos
      sxyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, charsize=csize
   endif

   if (keyword_set(sectarget)) then begin
      ypos = ypos + dypos
      sxyouts, xpos, ypos, 'SECTARGET = ' + sectarget, charsize=csize
   endif

end
;------------------------------------------------------------------------------

;+
; NAME:
;   plotspec
;
; PURPOSE:
;   Routine for plotting a single fiber from Princeton-1D spectro outputs
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
;   _EXTRA     - Kewords for SPLOT, such as XRANGE, YRANGE, THICK.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The data are read with READSPEC.  See the documentation for that
;   routine to see how to set environment variables that describe where
;   the data files are.
;
; EXAMPLES:
;   Plot the spectrum of plate 401, fiber #100 using the SPLOT plotting tool:
;     IDL> plotspec, 401, 100
;
;   The spectrum is shown in white, the errors in red (except masked points
;   are set to zero), and the best-fit eigenspectrum in blue. The mouse
;   buttons will zoom in (left), recenter (center), or zoom out (right).
;   The frame can be saved as a PostScript file by selecting File->WriteEPS
;   from the left-hand corner. 
;
;   Make the same plot, but boxcar-smooth the spectrum and limit the
;   wavelength range to [4000,5000] Angstroms:
;     IDL> plotspec, 401, 100, nsmooth=10, xrange=[5000,6000]
;
;   Some plates are observed on multiple nights. To select one of the two
;   observations of plate 306: 
;     IDL> plotspec, 306, 20, mjd=51690
;
; BUGS:
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
;-
;------------------------------------------------------------------------------
pro plotspec, plate, fiberid, mjd=mjd, nsmooth=nsmooth, $
 _EXTRA=KeywordsForSplot

   cspeed = 2.9979e5

   if (n_elements(plate) NE 1 OR n_elements(fiberid) NE 1) then $
    message, 'PLATE and FIBERID must be scalars'

   readspec, plate, fiberid, mjd=mjd, flux=objflux, flerr=objerr, $
    loglam=loglam, plug=plug, zans=zans
   wave = 10^loglam

   if (size(zans,/tname) EQ 'STRUCT') then begin
      synflux = synthspec(zans, loglam=loglam)
   endif else begin
      zans = 0
      synflux = 0
      print, 'WARNING: Redshift file not found'
   endelse

   if (keyword_set(nsmooth)) then begin
      objflux = smooth(objflux, nsmooth)
      if (keyword_set(synflux)) then $
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
    + '  MJD=' + strtrim(string(mjd),2)
   splot, wave, objflux, yrange=[ymin,ymax], $
    xtitle='Wavelength [Ang]', $
    ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
    title=title, charsize=csize, _EXTRA=KeywordsForSplot
   soplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
   soplot, wave, synflux, color='blue', _EXTRA=KeywordsForSplot, lw=2

   xpos = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
   dypos = 0.05 * (!y.crange[0] - !y.crange[1])
   ypos = !y.crange[1] + 1.5 * dypos

   if (keyword_set(zans)) then begin
      cz = zans.z * cspeed
      if (cz LT 1000) then $
        zstring = '  cz=' + string(cz,format='(f5.0)') + ' km/s' $
       else $
        zstring = '  z=' + string(zans.z,format='(f8.5)')
      sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
       charsize=csize

      ypos = ypos + dypos
      sxyouts, xpos, ypos, $
       TeXtoIDL('X^2_r =' + strtrim(string(zans.rchi2, format='(f6.2)'),2)), $
       charsize=csize
   endif

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

;+
; NAME:
;   plotspec
;
; PURPOSE:
;   Routine for plotting a single fiber from Princeton-1D spectro outputs
;
; CALLING SEQUENCE:
;   plotspec, plate, [ fiberid, mjd=, znum=, nsmooth=, psfile=, _EXTRA= ]
;
; INPUTS:
;   plate      - Plate number
;
; OPTIONAL INPUTS:
;   fiberid    - Fiber number(s); if not set, then plot all fibers for plate.
;   mjd        - MJD number; if not set, then select the most recent
;                data for this plate (largest MJD).
;   znum       - If set, then return not the best-fit redshift, but the
;                ZUM-th best-fit; e.g., set ZNUM=2 for second-best fit.
;   nsmooth    - If set, then boxcar smooth both the object and synthetic
;                spectra with a width equal to NSMOOTH.
;   psfile     - If set, then send plot to a PostScript file instead of
;                to the SPLOT interactive widget.  The PostScript file name
;                can be set explicitly, e.g. with PSFILE='test.ps'.  Or if
;                you simply set this as a flag, e.g. with /PSFILE, then the
;                default file name is spec-pppp-mmmmm-fff.ps,
;                where pppp=plate number, mmmmm=MJD, fff=fiber ID.
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
;   Loop through all the spectra for plate 401, interactively:
;     IDL> plotspec, 401
;
;   Plot all the spectra from plate 401 to a single PostScript file:
;     IDL> plotspec, 401, /psfile
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   readspec
;   soplot
;   splot
;   sdss_flagname()
;   sxyouts
;   synthspec()
;   textoidl()
;
; INTERNAL SUPPORT ROUTINES:
;   plotspec1
;
; REVISION HISTORY:
;   01-Sep-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro plotspec1, plate, fiberid, mjd=mjd, znum=znum, nsmooth=nsmooth, $
 psfile=psfile, xrange=passxr, yrange=passyr, noerase=noerase, $
 _EXTRA=KeywordsForSplot

   cspeed = 2.99792458e5

   readspec, plate, fiberid, mjd=mjd, znum=znum, flux=objflux, flerr=objerr, $
    loglam=loglam, plug=plug, zans=zans, synflux=synflux
   if (NOT keyword_set(objflux)) then begin
      print, 'Plate not found!!'
      return
   endif
   wave = 10^loglam

;   if (size(zans,/tname) EQ 'STRUCT') then begin
;      synflux = synthspec(zans, loglam=loglam)
;   endif else begin
;      zans = 0
;      synflux = 0
;      print, 'WARNING: Redshift file not found'
;   endelse

   if (keyword_set(nsmooth)) then begin
      if (nsmooth GT 1) then begin
         objflux = smooth(objflux, nsmooth)
         if (keyword_set(synflux)) then $
          synflux = smooth(synflux, nsmooth)
      endif
   endif

   primtarget = sdss_flagname('TARGET', plug.primtarget, /concat)
   sectarget = sdss_flagname('TTARGET', plug.sectarget, /concat)

   csize = 2
   textcolor = 'green'
   if (keyword_set(passyr)) then begin
      yrange = passyr
      ymin = yrange[0]
      ymax = yrange[1]
   endif else begin
      yrange = minmax(synflux)
      if (yrange[0] EQ yrange[1]) then yrange = minmax(objflux)
      ymin = (1.2 * yrange[0] - 0.2 * yrange[1]) < 0
      ymax = -0.2 * yrange[0] + 1.2 * yrange[1]
      if (ymax EQ ymin) then ymax = ymin + 1
      yrange = [ymin, ymax]
   endelse
   if (keyword_set(passxr)) then xrange = passxr

   title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(mjd),2)
   if (keyword_set(psfile)) then begin
      djs_plot, wave, objflux, xrange=xrange, yrange=yrange, $
       xtitle='Wavelength [Ang]', $
       ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
       title=title, charsize=csize, _EXTRA=KeywordsForSplot, /xstyle, /ystyle
      djs_oplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
      djs_oplot, wave, synflux, color='blue', _EXTRA=KeywordsForSplot, lw=2
   endif else begin
      if (NOT keyword_set(noerase)) then $
       splot, wave, objflux, xrange=xrange, yrange=yrange, $
        xtitle='Wavelength [Ang]', $
        ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
        title=title, charsize=csize, _EXTRA=KeywordsForSplot $
      else $
       soplot, wave, objflux, _EXTRA=KeywordsForSplot
      soplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
      soplot, wave, synflux, color='blue', _EXTRA=KeywordsForSplot, lw=2
   endelse

   xpos = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
   dypos = 0.05 * (!y.crange[0] - !y.crange[1])
   ypos = !y.crange[1] + 1.5 * dypos

   if (keyword_set(zans)) then begin
      cz = zans.z * cspeed
      if (abs(cz) LT 3000) then $
        zstring = '  cz=' + string(cz,format='(f6.0)') + ' km/s' $
       else $
        zstring = '  z=' + string(zans.z,format='(f8.5)')
      if (zans.zwarning NE 0) then $
       zstring = zstring + '  ZWARNING=' + strtrim(string(zans.zwarning),2) + ''
      if (keyword_set(znum)) then $
       zstring = zstring + ' (fit #' + strtrim(string(znum),2) + ')'

      if (keyword_set(psfile)) then $
       djs_xyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
        charsize=csize, color=textcolor $
      else $
       sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
        charsize=csize, color=textcolor

      ypos = ypos + dypos

      if (keyword_set(psfile)) then $
       djs_xyouts, xpos, ypos, $
        TeXtoIDL('X^2_r =' + strtrim(string(zans.rchi2, format='(f6.2)'),2)), $
        charsize=csize, color=textcolor $
      else $
       sxyouts, xpos, ypos, $
        TeXtoIDL('X^2_r =' + strtrim(string(zans.rchi2, format='(f6.2)'),2)), $
        charsize=csize, color=textcolor
   endif

   if (keyword_set(primtarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       djs_xyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=textcolor $
      else $
       sxyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=textcolor
   endif

   if (keyword_set(sectarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=textcolor $
      else $
       sxyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=textcolor
   endif

   return
end
;------------------------------------------------------------------------------
pro plotspec, plate, fiberid, mjd=mjd, znum=znum, nsmooth=nsmooth, $
 psfile=psfile, xrange=xrange, yrange=yrange, noerase=noerase, $
 _EXTRA=KeywordsForSplot

   if (n_params() LT 1) then begin
      print, 'Syntax - plotspec, plate, [ fiberid, mjd=, znum=, nsmooth=, $'
      print, '         psfile=, xrange=, yrange=, noerase=, $'
      print, '         _EXTRA=KeywordsForSplot'
      return
   endif

   if (n_elements(plate) NE 1) then $
    message, 'PLATE must be a scalar'
   if (NOT keyword_set(mjd)) then readspec, plate, mjd=mjd
   if (NOT keyword_set(mjd)) then begin
      print, 'Plate not found!!'
      return
   endif
   
   ;----------
   ; If writing to a PostScript file, then all plots are in the same file
   ; either if PSFILE is that file name, or if FIBERID is not specified
   ; (and then all 640 spectra are being plotted).

   if (NOT keyword_set(fiberid)) then begin
      fiberid = lindgen(640)+1L
      if (keyword_set(psfile)) then begin
         q_onefile = 1
         psfilename = string(plate, mjd, $
          format='("spec-",i4.4,"-",i5.5,".ps")')
      endif
   endif
   nfiber = n_elements(fiberid)
   if (size(psfile,/tname) EQ 'STRING' AND nfiber GT 1) then begin
      psfilename = psfile
      q_onefile = 1
   endif

   ;----------
   ; Loop over each plot

   ifiber = 0
   while (ifiber LT nfiber) do begin

      ;----------
      ; Open the PostScript file if appropriate

      if (keyword_set(psfile)) then begin
         if (NOT keyword_set(q_onefile)) then $
          psfilename = string(plate, mjd, fiberid[ifiber], $
           format='("spec-",i4.4,"-",i5.5,"-",i3.3,".ps")')

         if (NOT keyword_set(q_onefile) OR ifiber EQ 0) then begin
            dfpsplot, psfilename, /color, /square
         endif
      endif

      plotspec1, plate, fiberid[ifiber], mjd=mjd, znum=znum, $
       nsmooth=nsmooth, psfile=psfile, $
       xrange=xrange, yrange=yrange, noerase=noerase, _EXTRA=KeywordsForSplot

      if (keyword_set(psfile)) then begin
         if (NOT keyword_set(q_onefile) OR ifiber EQ nfiber-1) then dfpsclose
         ifiber = ifiber + 1
      endif else begin
         if (ifiber LT nfiber-1) then begin
            print, 'Press b=back one fiber'
            print, '      p=select new plate'
            print, '      f=select new fiber number'
            print, '      q=quit (and enter interactive mode for this plot)'
            print, '      s=change boxcar smoothing'
            print, '      x=change X plotting range'
            print, '      y=change Y plotting range'
            print, '      z=change which PCA-fit to overplot'
            print, '      any other key=forward'

            cc = strupcase(get_kbrd(1))
            case cc of
            'B': ifiber = (ifiber - 1) > 0
            'P': begin
                    read, plate, mjd, prompt='Enter new plate and MJD (enter 0 for unknown MJD): '
                    if (NOT keyword_set(mjd)) then readspec, plate, mjd=mjd
                    if (NOT keyword_set(mjd)) then begin
                       print, 'Plate not found!!'
                       return
                    endif
                 end
            'F': begin
                    read, newfiber, prompt='Enter new fiber number: '
                    nfiber = 640
                    fiberid = lindgen(nfiber) + 1L
                    ifiber = ((long(newfiber)-1) > 0) < (nfiber-1)
                 end
            'Q': ifiber = nfiber
            'S': begin
                    read, nsmooth, prompt='Enter boxcar smoothing width (0=none): '
                    nsmooth = long(nsmooth) > 0
                 end
            'X': begin
                    read, xmin, xmax, prompt='Enter new X range values (0 0=full range): '
                    if (xmin EQ 0 AND xmax EQ 0) then xrange = 0 $
                     else xrange = [xmin, xmax]
                 end
            'Y': begin
                    read, ymin, ymax, prompt='Enter new Y range values (0 0=full range): '
                    if (ymin EQ 0 AND ymax EQ 0) then yrange = 0 $
                     else yrange = [ymin, ymax]
                 end
            'Z': begin
                    read, znum, prompt='Enter 1=best redshift, 2=2nd best, ...: '
                    znum = long(znum) > 0
                 end
            else: ifiber = ifiber + 1
            endcase
         endif else begin
            ifiber = nfiber
         endelse
      endelse
   endwhile

   return
end
;------------------------------------------------------------------------------

;+
; NAME:
;   plotsignal
;
; PURPOSE:
;   Plot the signal (not S/N) seen on each plate+exposure+CCD.
;
; CALLING SEQUENCE:
;   plotsignal, plate, [ expnum, camname=, /addsky ]
;
; INPUTS:
;   plate      - Plate number
;
; OPTIONAL INPUTS:
;   expnum     - Exposure number(s); default to all exposures that have
;                an spFrame-b1 file on disk
;   camname    - Camera name(s); default to ['b2','r2','b1','r1']
;   addsky     - If set, then add the sky signal back in, rather than
;                using the sky-subtracted signal
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Most of our diagnostic plots show the S/N of the spectrographs.
;   These plots show the signal, and provide plots for each exposure.
;
;   One PostScript file is produced, with one page with one page per
;   exposure number.  The file names are:
;     spSignal-$PLATE.ps
;     spSignal-$PLATE-addsky.ps -- If /ADDSKY is specified
;
;   The spectro magnitudes are scaled to the 1st exposure on each
;   plate+MJD+camera combination.  There is a simple scaling for the
;   exposure time, but otherwise more-cloudy exposures will show as
;   negative (red) residuals.
;
;   This plots the signal after sky-subtraction, unless /ADDSKY is set.
;
; EXAMPLES:
;   Plot the signal seen on all exposures for plate 1489:
;     IDL> plotsignal, 1489
;
; BUGS:
;
; PROCEDURES CALLED:
;   bspline_valu()
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   filter_thru()
;   mrdfits()
;   splog
;   sxpar()
;   traceset2xy
;
; REVISION HISTORY:
;   30-Dec-2003  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro plotsignal, plate, expnum, camname=camname, addsky=addsky

   ;----------
   ; Check inputs and set defaults

   if (NOT keyword_set(plate)) then begin
      print, 'PLATE must be specified!'
      return
   endif
   if (NOT keyword_set(camname)) then camname = ['b2','r2','b1','r1']
   topdir = getenv('BOSS_SPECTRO_REDUX')
   csize = 1.1
   platestr = string(plate, format='(i4.4)')

   ;----------
   ; If EXPNUM is not set, then find all the exposure numbers that
   ; produced an spFrame file for the b1 CCD.

   if (NOT keyword_set(expnum)) then begin
      framefile = findfile(filepath('spFrame-b1-????????.fits*', $
       root_dir=topdir, subdir=platestr), count=ct)
      if (ct EQ 0) then begin
         print, 'No spFrame files found'
         return
      endif
      expnum = long(strmid(fileandpath(framefile),11,8))
   endif

   ; Open the plot file -- one page per exposure number
   if (keyword_set(addsky)) then $
    plotfile = string(plate, format='("spSignal-",i4.4,"-addsky.ps")') $
   else $
    plotfile = string(plate, format='("spSignal-",i4.4,".ps")')
   dfpsplot, plotfile, /color, /square

   last_mjd = lonarr(n_elements(camname))
   last_exptime = fltarr(n_elements(camname))
   magoffset = fltarr(n_elements(camname))

   for iexp=0, n_elements(expnum)-1 do begin
      splog, 'Working on exposure #', expnum[iexp]

      !p.multi = [0,2,2]
      for icam=0, n_elements(camname)-1 do begin
         speccolor = strmid(camname[icam],0,1)
         specid = fix(strmid(camname[icam],1,1))

         ; Read the spFrame file
         framefile = string(camname[icam], expnum[iexp], $
          format='("spFrame-",a2,"-",i8.8,".fits*")')
         framefile = (findfile(filepath(framefile, root_dir=topdir, $
          subdir=platestr)))[0]
         if (NOT keyword_set(framefile)) then message, 'File not found!!'
         fluximg = mrdfits(framefile, 0, hdr, /silent)
         fluxivar = mrdfits(framefile, 1, /silent)
         wset = mrdfits(framefile, 3, /silent)
         plug = mrdfits(framefile, 5, /silent)
         mjd = sxpar(hdr, 'MJD')
         exptime = sxpar(hdr, 'EXPTIME')

         ; Optionally add the sky flux back in
         if (keyword_set(addsky)) then $
          fluximg = fluximg + mrdfits(framefile, 6, /silent)

         ; Compute the wavelengths
         traceset2xy, wset, junk, loglam
         waveimg = 10d^loglam

         ; Read the spFluxcalib file -- if not found, then use a default
         fcalibfile = string(plate, mjd, camname[icam], $
          format='("spFluxcalib-",i4.4,"-",i5.5,"-",a2,".fits*")')
         fcalibfile = (findfile(filepath(fcalibfile, root_dir=topdir, $
          subdir=platestr)))[0]
         if (NOT keyword_set(fcalibfile)) then begin
            splog, 'Flux-calib file not found -- using defaults'
            fcalibfile = filepath('spFluxcalib-'+camname[icam]+'.fits', $
             root_dir=getenv('IDLSPEC2D_DIR'), subdir='examples')
         endif

         ; Compute the spectro-photometric calibration
         calibset = mrdfits(fcalibfile, 1, /silent)
         calibfac = bspline_valu(loglam, calibset)

         ; Integrate over the appropriate filter curve
         flambda2fnu = waveimg^2 / 2.99792e18
         specflux = filter_thru(fluximg/calibfac $
          * rebin(flambda2fnu,2048,320), $
          waveimg=waveimg, mask=(fluxivar LE 0))

         ifilt = (speccolor EQ 'b') ? 1 : 3 ; Either g-band or i-band
         specmag = -(48.6 - 2.5*17.) - 2.5*alog10(specflux[*,ifilt]>1e-15)
         photomag = plug.mag[ifilt]
         magdiff = photomag - specmag

         ; Scale the magnitudes to the 1st exposure on this plate+MJD+camera,
         ; using the exposure times.
         if (mjd NE last_mjd[icam]) then begin
            magoffset[icam] = median(magdiff)
            last_mjd[icam] = mjd
            last_exptime[icam] = exptime
         endif
         magdiff = magdiff - magoffset[icam] $
          - 2.5 * alog10(exptime / last_exptime[icam])

         symsize = abs(magdiff) * 1.
         qgood = photomag LT 20.5 AND specmag LT 21

         if (icam EQ 0) then $
          title=string(plate, mjd, strtrim(expnum[iexp],2), $
          format='("Plate ",i4, "  MJD ", i5, "  Exp# ", a)') $
         else $
          title = ''

         plot, [0], [0], /nodata, xchars=csize, ychars=csize, $
          xtitle='X [mm]', ytitle='Y [mm]', $
          xrange=[-320,320], yrange=[-320,320], xstyle=1, ystyle=1, $
          xmargin=xmargin, ymargin=ymargin, title=title
         ii = where(magdiff GT 0 AND qgood, ct)
         if (ct GT 1) then $
          djs_oplot, plug[ii].xfocal, plug[ii].yfocal, symsize=symsize[ii], $
           color='green', psym=2
         ii = where(magdiff LT 0 AND qgood, ct)
         if (ct GT 1) then $
          djs_oplot, plug[ii].xfocal, plug[ii].yfocal, symsize=symsize[ii], $
           color='red', psym=2
         xyouts, -300, 280, camname[icam], charsize=csize
         djs_oplot, [-300], [-280], psym=2, color='green'
         djs_oplot, [-300], [-300], psym=2, color='red'
         djs_xyouts, -300, -280, '  1 mag brighter', color='green'
         djs_xyouts, -300, -300, '  1 mag fainter', color='red'
      endfor
   endfor

   ; Close the plot file
   dfpsclose

   return
end
;------------------------------------------------------------------------------

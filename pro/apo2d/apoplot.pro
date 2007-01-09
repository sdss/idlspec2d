;+
; NAME:
;   apoplot
;
; PURPOSE:
;   Routine for plotting spectra from the Son-of-Spectro outputs at APO.
;
; CALLING SEQUENCE:
;   apoplot, plate, [ fiberid, mjd=, expnum=, nsmooth=, nmed=, psfile=, $
;    /magsort, /netimage, _EXTRA= ]
;
; INPUTS:
;   plate      - Plate number
;
; OPTIONAL INPUTS:
;   fiberid    - Fiber number(s); if not set, then plot all fibers for plate.
;   mjd        - MJD number; if not set, then select the most recent MJD
;                in the $SPECTROLOG_DIR directory.
;   expnum     - If set, then plot only these exposure numbers for this plate
;                rather than all exposure numbers for this plate.
;   nsmooth    - If set, then boxcar smooth the object spectra with a
;                width equal to NSMOOTH.
;   nmed       - If set, then median filter the object spectra with a
;                width equal to NMED.
;   psfile     - If set, then send plot to a PostScript file instead of
;                to the SPLOT interactive widget.  The PostScript file name
;                can be set explicitly, e.g. with PSFILE='test.ps'.  Or if
;                you simply set this as a flag, e.g. with /PSFILE, then the
;                default file name is spec-pppp-mmmmm-fff.ps,
;                where pppp=plate number, mmmmm=MJD, fff=fiber ID.
;   magsort    - If set and FIBERID is not, then plot all fibers from
;                the brightest object to the faintest.
;   netimage   - If set, then launch a Netscape browser with the object
;                image from Steve Kent's web site.  This only works if
;                Netscape is running and has permissions at the site
;                "http://sdssmosaic.fnal.gov:8015".
;                This is disabled if PSFILE is set.
;   _EXTRA     - Kewords for SPLOT, such as XRANGE, YRANGE, THICK.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The Son-of-Spectro outputs are first read from the file:
;     $SPECTROLOG_DIR/$MJD/logsheet-$MJD.fits
;   This then points to the other files that are read:
;     $SPECTROLOG_DIR/$MJD/fflat-$MJD-$PLATE-$EXPNUM-$CAMERA.fits
;     $SPECTROLOG_DIR/$MJD/wset-$MJD-$PLATE-$EXPNUM-$CAMERA.fits
;     $SPECTROLOG_DIR/$MJD/sci-$PLATE-$CAMERA-$EXPNUM.fits
;
;   If $SPECTROLOG_DIR is not set, then it is assumed to be
;     /data/spectro/spectrologs
;
;   The plotting range is set by the 5th and 95th percentiles of the data.
;   Note that there are a horrendous amount of cosmic rays in the data
;   as extracted by Son-of-Spectro.
;
;   The flux-calibration is very rudimentary, always using the same
;   four curves for the four cameras.
;
; EXAMPLES:
;   Plot the spectrum of plate 401, fiber #100 using the SPLOT plotting tool:
;     IDL> apoplot, 401, 100
;   The spectra from the first exposure (blue and red) are shown as white.
;   Other exposures are shown in other colors, and are labelled as such.
;   The mouse buttons will zoom in (left), recenter (center), or zoom out
;   (right).  The frame can be saved as a PostScript file by selecting
;   File->WriteEPS from the left-hand corner. 
;
;   Make the same plot, but boxcar-smooth the spectrum and limit the
;   wavelength range to [4000,5000] Angstroms:
;     IDL> apoplot, 401, 100, nsmooth=10, xrange=[5000,6000]
;
;   Some plates are observed on multiple nights. To select one of the two
;   observations of plate 306: 
;     IDL> apoplot, 306, 20, mjd=51690
;   This will only work if you have the Son-of-Spectro outputs for that
;   date on your disk.
;
;   Loop through all the spectra for plate 401, interactively:
;     IDL> apoplot, 401
;
;   Plot all the spectra from plate 401 to a single PostScript file:
;     IDL> apoplot, 401, /psfile
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/spFluxcalib-$CAMERA.fits
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   djs_median()
;   djs_oplot
;   djs_plot
;   fcalib_default()
;   get_mjd_dir()
;   soplot
;   splot
;   sdss_flagname()
;   sxyouts
;
; INTERNAL SUPPORT ROUTINES:
;   apoplot1
;
; REVISION HISTORY:
;   04-Dec-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apoplot1, plate, fiberid, mjd=mjd, expnum=allexpnum, nsmooth=nsmooth, $
 nmed=nmed, psfile=psfile, xrange=passxr, yrange=passyr, noerase=noerase, $
 netimage=netimage, _EXTRA=KeywordsForSplot

   common com_apoplot, mjddir, PPBIAS, PPFLAT, PPARC, PPSCIENCE

   ;----------
   ; Determine which spectrograph ID for this fiber

   if (fiberid LE 320) then begin
      specid = '1'
      column = [fiberid-1, fiberid-1]
   endif else begin
      specid = '2'
      column = [fiberid-321, fiberid-321]
   endelse

   ;----------
   ; Read the science spectra

   if (NOT keyword_set(PPSCIENCE)) then begin
      print, 'No science exposures for any plate on MJD ', mjd
      return
   endif
   sindx = where(PPSCIENCE.plate EQ plate $
    AND strmid(PPSCIENCE.camera,1,1) EQ specid, nscience)
   if (nscience EQ 0) then begin
      print, 'No science exposures for plate ', plate
      return
   endif
   for iscience=0, nscience-1 do begin
      thisfile = filepath(PPSCIENCE[sindx[iscience]].scifile, root_dir=mjddir)
      if (iscience EQ 0) then begin
         objsub = mrdfits(thisfile, 0, range=column, /silent)
         objsubivar = mrdfits(thisfile, 1, range=column, /silent)
      endif else begin
         objsub = [[objsub], [mrdfits(thisfile, 0, range=column, /silent)]]
         objsubivar = [[objsubivar], [mrdfits(thisfile, range=column, /silent)]]
      endelse
   endfor

   allcams = PPSCIENCE[sindx].camera
   allcams = allcams[ uniq(allcams, sort(allcams)) ]
   cindx = lonarr(nscience)
   for icam=0, n_elements(allcams)-1 do begin
      ; Select which science frames correspond to this plate+camera
      jscience = where(PPSCIENCE[sindx].camera EQ allcams[icam])
      cindx[jscience] = icam

      ; Read the flat for this plate+camera
;      j = where(PPFLAT.plate EQ plate AND PPFLAT.camera EQ allcams[icam])
;      j = (reverse(j))[0] ; Select the last one
;      thisfile = filepath(PPFLAT[j].tsetfile, root_dir=mjddir)
;      flat = mrdfits(thisfile, 0, range=column, /silent)
;      ; Apply this flat to corresponding science frames
;      for jj=0, n_elements(jscience)-1 do begin
;         objsub[*,jscience[jj]] = objsub[*,jscience[jj]] / flat ; zeros?
;         objsubivar[*,jscience[jj]] = objsubivar[*,jscience[jj]] * flat^2
;      endfor

      ; Only need to read the plugmap once
      j = where(PPFLAT.plate EQ plate AND PPFLAT.camera EQ allcams[icam])
      j = (reverse(j))[0] ; Select the last one
      if (icam EQ 0) then $
       plug = mrdfits(filepath(PPFLAT[j].tsetfile, root_dir=mjddir), $
        3, range=column, /silent)

      ; Read the arc for this plate+camera
      j = where(PPARC.plate EQ plate AND PPARC.camera EQ allcams[icam])
      j = (reverse(j))[0] ; Select the last one
      thisfile = filepath(PPARC[j].wsetfile, root_dir=mjddir)
      wset = mrdfits(thisfile, 1, /silent)
      traceset2xy, wset, xarc, yarc
      thisloglam = yarc[*,column[0]]
      if (icam EQ 0) then loglam = thisloglam $
       else loglam = [[loglam], [thisloglam]]

      ; Read and apply a canonical flux-calibration vector
      calibfac = fcalib_default(allcams[icam], yarc[*,column[0]], 1.)
      qgood = calibfac GT 0
      for jj=0, n_elements(jscience)-1 do begin
         objsub[*,jscience[jj]] = qgood * objsub[*,jscience[jj]] $
          / (PPSCIENCE[jscience[jj]].exptime * calibfac + (qgood EQ 0))
      endfor

      ; Mask out any regions where we don't know the flux-calibration vector
      for jj=0, n_elements(jscience)-1 do $
       objsubivar[*,jscience[jj]] = objsubivar[*,jscience[jj]] * (calibfac NE 0)

   endfor

   wave = 10^loglam

   ;----------
   ; Choose plot colors based upon exposure numbers

   if (NOT keyword_set(allexpnum)) then begin
      allexpnum = PPSCIENCE[sindx].expnum
      allexpnum = allexpnum[ uniq(allexpnum, sort(allexpnum)) ]
   endif
   nexp = n_elements(allexpnum)

   ;----------
   ; Interpolate over bad pixels and smooth if specified

   objsub = djs_maskinterp(objsub, objsubivar EQ 0, iaxis=0, /const)
   for iscience=0, nscience-1 do begin
      if (keyword_set(nmed)) then begin
         if (nmed GT 1) then begin
            objsub[*,iscience] = djs_median(objsub[*,iscience], width=nmed, $
             boundary='reflect')
         endif
      endif
      if (keyword_set(nsmooth)) then begin
         if (nsmooth GT 1) then begin
            objsub[*,iscience] = smooth(objsub[*,iscience], nsmooth)
         endif
      endif
   endfor

   ;----------
   ; Make the actual plots

   csize = 1.75
   textcolor = 'default'
   colorvec = ['default', 'red', 'green', 'blue', 'magenta', 'cyan']
   title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(mjd),2)
   primtarget = sdss_flagname('TARGET', plug.primtarget, /concat)
   sectarget = sdss_flagname('TTARGET', plug.sectarget, /concat)

   if (keyword_set(passxr)) then xrange = passxr $
    else xrange = minmax(wave)
   if (keyword_set(passyr)) then begin
      yrange = passyr
   endif else begin
;      yrange = minmax(objsub)
      ; Use the 5th and 95th percentiles for setting the Y plot range
      isort = sort(objsub)
      ymin = objsub[isort[long(0.05*n_elements(objsub))]]
      ymax = objsub[isort[long(0.95*n_elements(objsub))]]
      yrange = [ymin, ymax]
      ymin = (1.2 * yrange[0] - 0.2 * yrange[1]) < 0
      ymax = -0.2 * yrange[0] + 1.2 * yrange[1]
      if (ymax EQ ymin) then ymax = ymin + 1
      yrange = [ymin, ymax]
   endelse

   if (keyword_set(psfile)) then $
    djs_plot, [0], [0], /nodata, charsize=csize, $
     xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
     xtitle='Wavelength [Ang]', ytitle='Flux [electrons]', title=title, $
     _EXTRA=KeywordsForSplot $
   else $
    splot, [0], [0], /nodata, charsize=csize, $
     xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
     xtitle='Wavelength [Ang]', ytitle='Flux [electrons]', title=title, $
     _EXTRA=KeywordsForSplot
   for iscience=0, nscience-1 do begin
      ; Only plot this exposure if it is one listed in ALLEXPNUM.
      ithis = (where(PPSCIENCE[sindx[iscience]].expnum EQ allexpnum))[0]
      if (ithis NE -1) then begin
         thiscolor = colorvec[ithis MOD n_elements(colorvec)]
         if (keyword_set(psfile)) then $
          djs_oplot, wave[*,cindx[iscience]], objsub[*,iscience], $
           color=thiscolor, _EXTRA=KeywordsForSplot $
         else $
          soplot, wave[*,cindx[iscience]], objsub[*,iscience], $
           color=thiscolor, _EXTRA=KeywordsForSplot
      endif
   endfor

   xpos = 0.80 * !x.crange[0] + 0.20 * !x.crange[1]
   dypos = 0.05 * (!y.crange[1] - !y.crange[0])
   ypos = !y.crange[0] + 0.0 * dypos

   if (keyword_set(primtarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=djs_icolor(textcolor) $
      else $
       sxyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=textcolor
   endif

   if (keyword_set(sectarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=djs_icolor(textcolor) $
      else $
       sxyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=textcolor
   endif

   for iexp=0, nexp-1 do begin
      ypos = ypos + dypos
      thiscolor = djs_icolor(colorvec[iexp MOD n_elements(colorvec)])
      if (keyword_set(psfile)) then $
       djs_xyouts, xpos, ypos, $
       'Exposure #' + strtrim(string(allexpnum[iexp]),2), $
       charsize=2.0, color=thiscolor $
      else $
       sxyouts, xpos, ypos, $
        'Exposure #' + strtrim(string(allexpnum[iexp]),2), $
        charsize=2.0, color=thiscolor
   endfor

   if (keyword_set(netimage) AND NOT keyword_set(psfile)) then begin
      netstring = 'http://sdssmosaic.fnal.gov:8015/template/tsSingle.tml?run=' $
       + strtrim(string(plug.objid[0]),2) $
       + '&camcol=' + strtrim(string(plug.objid[2]),2) $
       + '&field=' + strtrim(string(plug.objid[3]),2) $
       + '&ra=' + strtrim(string(plug.ra),2) $
       + '&dec=' + strtrim(string(plug.dec),2)
      spawn, '\netscape -remote "openURL(' + netstring + ')"'
   endif

   return
end

;------------------------------------------------------------------------------
pro apoplot, plate, fiberid, mjd=mjd, expnum=expnum, nsmooth=nsmooth, $
 nmed=nmed, psfile=psfile, noerase=noerase, xrange=xrange, yrange=yrange, $
 magsort=magsort, _EXTRA=KeywordsForSplot

   common com_apoplot, mjddir, PPBIAS, PPFLAT, PPARC, PPSCIENCE

   if (n_params() LT 1) then begin
      print, 'Syntax - apoplot, plate, [ fiberid, mjd=, nsmooth=, $'
      print, '         nmed=, psfile=, xrange=, yrange=, /noerase, $'
      print, '         /netimage, _EXTRA=KeywordsForSplot'
      return
   endif

   quiet = !quiet
   !quiet = 1

   if (n_elements(plate) NE 1) then $
    message, 'PLATE must be a scalar'

   ;----------
   ; If MJD is not specified, then find the most recent MJD for output files

   spectrolog_dir = getenv('SPECTROLOG_DIR')
   if (NOT keyword_set(spectrolog_dir)) then $
    spectrolog_dir = '/data/spectro/spectrologs'

   if (NOT keyword_set(mjd)) then begin
      mjdlist = get_mjd_dir(spectrolog_dir, mjstart=1, mjend=99999, mjd='?????')
      mjd = (reverse(mjdlist[sort(mjdlist)]))[0]
      splog, 'Selecting MJD=', mjd, ' (override this with MJD keyword)'
   endif

   ;----------
   ; Read the log file

   mjdstr = string(mjd, format='(i5.5)')
   mjddir = concat_dir(spectrolog_dir, mjdstr)

   logfile = 'logfile-' + mjdstr + '.fits'
   logfile = filepath(logfile, root_dir=mjddir)
   if (NOT keyword_set(findfile(logfile))) then begin
      print, 'Unable to find logfile '+logfile
      !quiet = quiet
      return
   endif

   ; Read in all the HDU's in the log file as structures
   PPBIAS = mrdfits(logfile, 1, /silent)
   PPFLAT = mrdfits(logfile, 2, /silent)
   PPARC = mrdfits(logfile, 3, /silent)
   PPSCIENCE = mrdfits(logfile, 4, /silent)

   ;----------
   ; If writing to a PostScript file, then all plots are in the same file
   ; either if PSFILE is that file name, or if FIBERID is not specified
   ; (and then all 640 spectra are being plotted).

   if (NOT keyword_set(fiberid)) then begin
      if (keyword_set(magsort) AND keyword_set(PPSCIENCE)) then begin
         sindx1 = (where(PPSCIENCE.plate EQ plate $
          AND (strmid(PPSCIENCE.camera,1,1) EQ '1')))[0]
         sindx2 = (where(PPSCIENCE.plate EQ plate $
          AND (strmid(PPSCIENCE.camera,1,1) EQ '2')))[0]
         if (sindx2 EQ -1) then begin
            fiberid = lindgen(320)+1
            fibermag = PPSCIENCE[sindx1].fibermag
         endif else if (sindx1 EQ -1) then begin
            fiberid = lindgen(320)+321
            fibermag = PPSCIENCE[sindx2].fibermag
         endif else begin
            fiberid = lindgen(640)
            fibermag = [PPSCIENCE[sindx1].fibermag, PPSCIENCE[sindx2].fibermag]
         endelse
         fiberid = fiberid[ sort(fibermag) ]
      endif
      if (NOT keyword_set(fiberid)) then $
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

      apoplot1, plate, fiberid[ifiber], mjd=mjd, expnum=expnum, $
       nmed=nmed, nsmooth=nsmooth, psfile=psfile, $
       xrange=xrange, yrange=yrange, noerase=noerase, netimage=netimage, $
       _EXTRA=KeywordsForSplot

      if (keyword_set(psfile)) then begin
         if (NOT keyword_set(q_onefile) OR ifiber EQ nfiber-1) then dfpsclose
         ifiber = ifiber + 1
      endif else begin
         if (ifiber LT nfiber-1) then begin
            if (keyword_set(nsmooth)) then $
             sstring = ' (currently=' + strtrim(string(nsmooth),2) + ')' $
            else $
             sstring = ''

            if (keyword_set(nmed)) then $
             mstring = ' (currently=' + strtrim(string(nmed),2) + ')' $
            else $
             mstring = ''

            print, 'Press b=back one fiber'
            print, '      p=select new plate'
            print, '      f=select new fiber number'
            print, '      q=quit (and enter interactive mode for this plot)'
            print, '      m=change median filtering' + mstring
            print, '      s=change boxcar smoothing' + sstring
            print, '      x=change X plotting range'
            print, '      y=change Y plotting range'
            print, '      any other key=forward'

            cc = strupcase(get_kbrd(1))
            case cc of
            'B': ifiber = (ifiber - 1) > 0
            'P': begin
                    read, plate, prompt='Enter new plate (same MJD assumed): '
                 end
            'F': begin
                    read, newfiber, prompt='Enter new fiber number: '
                    nfiber = 640
                    fiberid = lindgen(nfiber) + 1L
                    ifiber = ((long(newfiber)-1) > 0) < (nfiber-1)
                 end
            'Q': ifiber = nfiber
            'M': begin
                    read, nmed, prompt='Enter median filtering width (0=none): '
                    nmed = long(nmed) > 0
                 end
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
            else: ifiber = ifiber + 1
            endcase
         endif else begin
            ifiber = nfiber
         endelse
      endelse
   endwhile

   !quiet = quiet
   return
end
;------------------------------------------------------------------------------
